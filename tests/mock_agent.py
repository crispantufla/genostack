"""Minimal ACP agent used in tests: answers each prompt with a JSON block matching genostack's schemas.

Run manually:  python tests/mock_agent.py   (speaks ACP over stdio)
"""
from __future__ import annotations

import asyncio
import json
import re
import sys
import uuid

from acp import PROTOCOL_VERSION, run_agent
from acp.helpers import session_notification, update_agent_message_text
from acp.schema import (AgentCapabilities, Implementation, InitializeResponse, NewSessionResponse, PromptResponse)


class MockAgent:
    def __init__(self):
        self.conn = None

    def on_connect(self, conn):
        self.conn = conn

    async def initialize(self, protocol_version, client_capabilities=None, client_info=None, **kw):
        return InitializeResponse(protocolVersion=PROTOCOL_VERSION, agentCapabilities=AgentCapabilities(),
                                  agentInfo=Implementation(name="mock-agent", version="0.0.1"))

    async def authenticate(self, method_id, **kw):
        return None

    async def new_session(self, cwd, mcp_servers=None, **kw):
        return NewSessionResponse(sessionId=str(uuid.uuid4()))

    async def load_session(self, *a, **kw):
        return None

    async def set_session_mode(self, *a, **kw):
        return None

    async def cancel(self, session_id, **kw):
        return None

    async def prompt(self, session_id, prompt, **kw):
        text = "".join(getattr(b, "text", "") for b in prompt)
        m = re.search(r'"module":\s*"([^"]+)"', text)
        if "core_stack" in text:  # synthesis prompt
            payload = {
                "executive_summary": "Resumen de prueba.", "genetic_profile_summary": "Perfil de prueba.",
                "core_stack": [{"name": "Creatina", "why": "test", "dose": "5 g", "timing": "diario", "evidence_grade": "A", "priority": 1, "category": "suplemento", "genes": ["MTHFR"]}],
                "advanced_experimental": [], "avoid": ["nada"], "conflicts": ["ninguno"], "introduction_order": ["fase 1"],
                "labs": ["homocisteína"], "confirm_with_sequencing": [], "unexpected_findings": [], "caveats": "mock",
            }
        else:
            mod = m.group(1) if m else "unknown"
            payload = {
                "module": mod, "interpretation": f"Interpretación mock para {mod}.",
                "interventions": [{"name": "Riboflavina", "name_en": "riboflavin", "category": "suplemento", "target_genes": ["MTHFR"], "rsids": ["rs1801133"],
                                   "mechanism": "cofactor FAD", "expected_effect": "↓ homocisteína", "effect_magnitude": "−40 % en TT",
                                   "evidence_grade": "A", "key_studies": ["McNulty 2006 Circulation (PMID 16380544)"], "dose": "1.6 mg/d",
                                   "risks_interactions": "ninguno", "novelty": "establecido", "legal_status": "OTC", "priority": 1}],
                "avoid": [], "labs": ["homocisteína"], "open_questions": [],
            }
        chunk = "Aquí tienes:\n```json\n" + json.dumps(payload, ensure_ascii=False) + "\n```"
        for i in range(0, len(chunk), 400):
            await self.conn.session_update(session_id=session_id, update=update_agent_message_text(chunk[i:i + 400]))
        return PromptResponse(stopReason="end_turn")

    async def ext_method(self, method, params):
        return {}

    async def ext_notification(self, method, params):
        return None


if __name__ == "__main__":
    asyncio.run(run_agent(MockAgent()))
