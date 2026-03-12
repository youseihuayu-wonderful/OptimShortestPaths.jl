"""
ChemPath Agent Engine — Claude API conversation loop with tool use.

Manages multi-turn conversations, tool dispatch, and response formatting.
"""

from __future__ import annotations

import os
from typing import Generator

from chempath.agent.tools import TOOL_DEFINITIONS, ChemPathToolExecutor

SYSTEM_PROMPT = """\
You are ChemPath, an expert drug screening assistant. You help researchers find \
and evaluate drug candidates for protein targets using graph-based optimization.

You have access to a database of real bioactivity data from ChEMBL containing \
thousands of compounds tested against anti-cancer targets (EGFR, BRAF, ALK, HER2, \
VEGFR2, ABL1, MET, FGFR1, PI3Ka, mTOR).

Your capabilities:
- Screen and rank compounds for any target using shortest-path graph optimization
- Compare strategies: balanced (efficacy + safety), efficacy-first, safety-first
- Compute Pareto-optimal trade-offs between efficacy and toxicity
- Run sensitivity analysis to test robustness of recommendations
- Look up detailed info on specific compounds

When presenting results:
- Always explain what the numbers mean in plain language
- Highlight the confidence level (HIGH/MEDIUM/LOW) and explain why
- Note limitations: these are computational rankings, not clinical recommendations
- If toxicity data is unavailable, mention that safety analysis is limited

Be concise but thorough. Use the tools to get data before answering questions \
about specific compounds or targets.\
"""

DEFAULT_MODEL = "claude-sonnet-4-20250514"


class ChemPathAgent:
    """Multi-turn conversation agent with Claude tool use."""

    def __init__(
        self,
        api_key: str | None = None,
        model: str = DEFAULT_MODEL,
        data_path: str | None = None,
    ):
        try:
            import anthropic
        except ImportError:
            raise ImportError(
                "anthropic package required. Install with: uv add anthropic"
            )

        self.client = anthropic.Anthropic(api_key=api_key or os.environ.get("ANTHROPIC_API_KEY"))
        self.model = model
        self.executor = ChemPathToolExecutor(data_path) if data_path else ChemPathToolExecutor()
        self.messages: list[dict] = []

    def chat(self, user_message: str) -> str:
        """Send a message and get a response, handling tool calls automatically."""
        self.messages.append({"role": "user", "content": user_message})

        # Agentic loop: keep calling Claude until we get a final text response
        while True:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=4096,
                system=SYSTEM_PROMPT,
                tools=TOOL_DEFINITIONS,
                messages=self.messages,
            )

            # Collect all content blocks
            assistant_content = response.content
            self.messages.append({"role": "assistant", "content": assistant_content})

            # If stop_reason is "end_turn", extract text and return
            if response.stop_reason == "end_turn":
                text_parts = [
                    block.text for block in assistant_content
                    if hasattr(block, "text")
                ]
                return "\n".join(text_parts) if text_parts else "(No response)"

            # If stop_reason is "tool_use", execute tools and continue
            if response.stop_reason == "tool_use":
                tool_results = []
                for block in assistant_content:
                    if block.type == "tool_use":
                        result = self.executor.execute(block.name, block.input)
                        tool_results.append({
                            "type": "tool_result",
                            "tool_use_id": block.id,
                            "content": result,
                        })
                self.messages.append({"role": "user", "content": tool_results})
                continue

            # Unexpected stop reason
            text_parts = [
                block.text for block in assistant_content
                if hasattr(block, "text")
            ]
            return "\n".join(text_parts) if text_parts else "(Unexpected stop)"

    def chat_stream(self, user_message: str) -> Generator[str, None, None]:
        """Stream a response, yielding text chunks as they arrive."""
        self.messages.append({"role": "user", "content": user_message})

        while True:
            collected_content = []
            current_tool_uses = {}

            with self.client.messages.stream(
                model=self.model,
                max_tokens=4096,
                system=SYSTEM_PROMPT,
                tools=TOOL_DEFINITIONS,
                messages=self.messages,
            ) as stream:
                for event in stream:
                    if hasattr(event, "type"):
                        if event.type == "content_block_start":
                            block = event.content_block
                            if hasattr(block, "text"):
                                collected_content.append(block)
                            elif block.type == "tool_use":
                                current_tool_uses[block.id] = {
                                    "id": block.id,
                                    "name": block.name,
                                    "input_json": "",
                                }
                        elif event.type == "content_block_delta":
                            delta = event.delta
                            if hasattr(delta, "text"):
                                yield delta.text
                            elif hasattr(delta, "partial_json"):
                                # Accumulate tool input JSON
                                for tu in current_tool_uses.values():
                                    tu["input_json"] += delta.partial_json

                response = stream.get_final_message()

            self.messages.append({"role": "assistant", "content": response.content})

            if response.stop_reason == "end_turn":
                return

            if response.stop_reason == "tool_use":
                tool_results = []
                for block in response.content:
                    if block.type == "tool_use":
                        yield f"\n[Calling {block.name}...]\n"
                        result = self.executor.execute(block.name, block.input)
                        tool_results.append({
                            "type": "tool_result",
                            "tool_use_id": block.id,
                            "content": result,
                        })
                self.messages.append({"role": "user", "content": tool_results})
                continue

            return

    def reset(self):
        """Clear conversation history."""
        self.messages = []
