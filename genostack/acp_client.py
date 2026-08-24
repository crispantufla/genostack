"""Thin ACP (Agent Client Protocol) client: spawn an agent, run research prompts, collect text.

We are the *client*; the agent (Claude Code via claude-agent-acp, Gemini CLI, codex-acp, ...) is a
subprocess speaking JSON-RPC over stdio. Permissions: read/search/fetch/think allowed, anything that
writes or executes is rejected. Each prompt runs in its own session; several sessions can run concurrently.
"""
from __future__ import annotations

import asyncio
import json
import os
import re
import shlex
import shutil
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Callable

from acp import PROTOCOL_VERSION, spawn_agent_process, text_block
from acp.interfaces import Client
from acp.exceptions import RequestError
from acp.schema import AllowedOutcome, DeniedOutcome, RequestPermissionResponse


class AgentAuthError(RuntimeError):
    """The agent cannot authenticate; there is no point retrying other tasks."""


class AgentQuotaError(RuntimeError):
    """The account hit its usage limit; retrying immediately is pointless."""


# transient server-side failures worth retrying with backoff
TRANSIENT_RE = re.compile(r"\b(429|500|502|503|504|529)\b|overloaded|rate.?limit|temporar|timed? ?out|"
                          r"connection (reset|closed|error)|stopped arriving|internal server|sin texto", re.I)
QUOTA_RE = re.compile(r"usage limit|quota (exceeded|reached)|out of credit|insufficient (credit|balance)|"
                      r"limit reached|upgrade to continue", re.I)

ALLOWED_KINDS = {"read", "search", "fetch", "think", "other"}
DENIED_TITLE_RE = re.compile(r"\b(bash|shell|powershell|write|edit|multiedit|notebookedit|rm\b|delete|move|execute|terminal)\b", re.I)

ProgressCb = Callable[[str, str], None]   # (label, status-text)


class _Client(Client):
    def __init__(self, owner: "AgentRunner"):
        self.owner = owner

    async def request_permission(self, session_id, tool_call, options, **kwargs):
        kind = getattr(tool_call, "kind", None) or "other"
        title = getattr(tool_call, "title", "") or ""
        allow = kind in ALLOWED_KINDS and not DENIED_TITLE_RE.search(title)
        if kind == "other" and DENIED_TITLE_RE.search(title):
            allow = False
        want = ("allow_once", "allow_always") if allow else ("reject_once", "reject_always")
        for w in want:
            for o in options:
                if o.kind == w:
                    if allow:
                        return RequestPermissionResponse(outcome=AllowedOutcome(outcome="selected", optionId=o.option_id))
                    return RequestPermissionResponse(outcome=AllowedOutcome(outcome="selected", optionId=o.option_id))
        return RequestPermissionResponse(outcome=DeniedOutcome(outcome="cancelled"))

    async def session_update(self, session_id, update, **kwargs):
        st = self.owner.sessions.get(session_id)
        if st is None:
            return
        su = getattr(update, "session_update", "")
        if su == "agent_message_chunk":
            c = getattr(update, "content", None)
            txt = getattr(c, "text", None)
            if txt:
                st["text"].append(txt)
                st["chars"] += len(txt)
        elif su in ("tool_call", "tool_call_update"):
            title = getattr(update, "title", None)
            if title:
                st["last_tool"] = title
                st["tools"] += 1 if su == "tool_call" else 0
        elif su == "agent_thought_chunk":
            st["thoughts"] += 1
        self.owner._progress(st)

    async def write_text_file(self, *a, **k):
        raise RuntimeError("write denied")

    async def read_text_file(self, *a, **k):
        raise RuntimeError("read denied")

    async def create_terminal(self, *a, **k):
        raise RuntimeError("terminal denied")

    async def terminal_output(self, *a, **k):
        raise RuntimeError("terminal denied")

    async def release_terminal(self, *a, **k):
        return None

    async def wait_for_terminal_exit(self, *a, **k):
        raise RuntimeError("terminal denied")

    async def kill_terminal(self, *a, **k):
        return None

    async def create_elicitation(self, *a, **k):
        raise RuntimeError("elicitation unsupported")

    async def complete_elicitation(self, *a, **k):
        return None

    async def ext_method(self, method, params):
        return {}

    async def ext_notification(self, method, params):
        return None


class AgentRunner:
    """Usage:
        async with AgentRunner("npx -y @agentclientprotocol/claude-agent-acp", progress=cb) as ar:
            text = await ar.ask("...", label="methylation")
    """

    def __init__(self, command: str, progress: ProgressCb | None = None, concurrency: int = 3, log_dir: Path | None = None):
        self.command = command
        self.progress_cb = progress
        self.sem = asyncio.Semaphore(max(1, concurrency))
        self.sessions: dict[str, dict] = {}
        self._cm = None
        self.conn = None
        self.proc = None
        self.log_dir = log_dir or Path(tempfile.mkdtemp(prefix="genostack-agent-"))
        self._stderr_f = None
        self.agent_info: str = ""

    async def __aenter__(self):
        parts = [p.strip('"') for p in shlex.split(self.command, posix=(os.name != "nt"))]
        exe = shutil.which(parts[0]) or parts[0]
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self._stderr_f = open(self.log_dir / "agent_stderr.log", "ab")
        # pass through auth/proxy variables (the SDK trims the environment to a safe minimum by default)
        env = {k: v for k, v in os.environ.items()
               if k.startswith(("ANTHROPIC_", "CLAUDE_CODE_OAUTH", "GEMINI_", "GOOGLE_", "OPENAI_", "HTTP_PROXY", "HTTPS_PROXY", "NO_PROXY", "http_proxy", "https_proxy"))
               or k in ("HOME", "COMSPEC", "ProgramFiles", "ProgramData", "NVM_HOME", "NVM_SYMLINK", "NODE_PATH", "npm_config_prefix")}
        self._cm = spawn_agent_process(_Client(self), exe, *parts[1:], cwd=str(self.log_dir), env=env,
                                       transport_kwargs={"stderr": self._stderr_f.fileno(), "limit": 32 * 1024 * 1024})
        self.conn, self.proc = await self._cm.__aenter__()
        try:
            init = await asyncio.wait_for(self.conn.initialize(protocol_version=PROTOCOL_VERSION), timeout=120)
        except Exception as e:  # noqa: BLE001
            self._stderr_f.flush()
            tail = (self.log_dir / "agent_stderr.log").read_text(encoding="utf-8", errors="replace")[-600:]
            raise RuntimeError(f"No se pudo inicializar el agente ACP `{self.command}`: {e}. stderr: {tail.strip()}") from e
        info = getattr(init, "agent_info", None)
        self.agent_info = f"{getattr(info, 'name', '')} {getattr(info, 'version', '')}".strip() or "agente ACP"
        auth = getattr(init, "auth_methods", None) or []
        self._auth_methods = [getattr(a, "id", "") for a in auth]
        return self

    async def __aexit__(self, *exc):
        try:
            if self._cm:
                await self._cm.__aexit__(*exc)
        finally:
            if self._stderr_f:
                self._stderr_f.close()

    def _progress(self, st: dict):
        if self.progress_cb:
            self.progress_cb(st["label"], f"{st['chars']//1000}k chars · {st['tools']} tools · {st['last_tool'][:50]}")

    async def ask(self, prompt: str, label: str = "task", timeout: float = 1200.0, attempts: int = 4) -> str:
        """Run one prompt in its own session, retrying transient server failures with backoff."""
        last: Exception | None = None
        for attempt in range(attempts):
            async with self.sem:
                try:
                    return await self._ask_once(prompt, label if attempt == 0 else f"{label}.try{attempt + 1}", timeout)
                except (AgentAuthError, AgentQuotaError):
                    raise
                except RuntimeError as e:
                    if not TRANSIENT_RE.search(str(e)) or attempt == attempts - 1:
                        raise
                    last = e
            wait = min(120, 15 * 2 ** attempt)
            if self.progress_cb:
                self.progress_cb(label, f"error transitorio, reintento en {wait}s ({attempt + 1}/{attempts - 1})")
            await asyncio.sleep(wait)
        raise last if last else RuntimeError(f"sesión {label} sin resultado")

    async def _ask_once(self, prompt: str, label: str, timeout: float) -> str:
        sess = await asyncio.wait_for(self.conn.new_session(cwd=str(self.log_dir), mcp_servers=[]), timeout=120)
        sid = sess.session_id
        st = {"label": label, "text": [], "chars": 0, "tools": 0, "last_tool": "", "thoughts": 0, "t0": time.time()}
        self.sessions[sid] = st
        self._progress(st)
        try:
            resp = await asyncio.wait_for(self.conn.prompt(session_id=sid, prompt=[text_block(prompt)]), timeout=timeout)
            stop = getattr(resp, "stop_reason", "")
        except RequestError as e:
            msg = str(e)
            if "authenticate" in msg.lower() or "oauth" in msg.lower() or "api key" in msg.lower():
                raise AgentAuthError(
                    f"El agente ACP no pudo autenticarse ({msg[:160]}). Inicia sesión en Claude Code desde una terminal "
                    "(`claude` → /login) o exporta ANTHROPIC_API_KEY antes de ejecutar genostack.") from e
            if QUOTA_RE.search(msg):
                raise AgentQuotaError(
                    f"Se agotó el límite de uso de la cuenta ({msg[:160]}). Espera a que se restablezca (o usa "
                    "ANTHROPIC_API_KEY) y relanza con --resume para no repetir los módulos ya investigados.") from e
            raise RuntimeError(f"Error del agente en la sesión {label}: {msg[:300]}") from e
        except asyncio.TimeoutError:
            stop = "timeout"
            try:
                await self.conn.cancel(session_id=sid)
            except Exception:  # noqa: BLE001
                pass
        text = "".join(st["text"])
        (self.log_dir / f"{label}.prompt.md").write_text(prompt, encoding="utf-8")
        (self.log_dir / f"{label}.response.md").write_text(text, encoding="utf-8")
        if not text.strip():
            raise RuntimeError(f"El agente terminó con stop_reason={stop} sin texto (sesión {label}) "
                               "— suele indicar corte del servidor o límite de uso; se reintentará")
        st["stop"] = stop
        return text


def extract_json(text: str) -> Any:
    """Pull the first JSON object out of an agent reply (fenced or bare)."""
    m = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", text, re.S)
    candidates = [m.group(1)] if m else []
    # bare: first '{' to last '}'
    i, j = text.find("{"), text.rfind("}")
    if i != -1 and j > i:
        candidates.append(text[i:j + 1])
    for c in candidates:
        try:
            return json.loads(c)
        except json.JSONDecodeError:
            # try to fix trailing commas
            try:
                return json.loads(re.sub(r",(\s*[}\]])", r"\1", c))
            except json.JSONDecodeError:
                continue
    raise ValueError("No se encontró JSON válido en la respuesta del agente")
