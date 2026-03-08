"""
ChemPath Chainlit Chat UI.

Run with: uv run chainlit run chempath/ui/app.py
"""

import chainlit as cl

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

MODEL = "claude-sonnet-4-20250514"

executor = ChemPathToolExecutor()


@cl.on_chat_start
async def on_start():
    """Initialize the chat session."""
    import anthropic

    client = anthropic.AsyncAnthropic()
    cl.user_session.set("client", client)
    cl.user_session.set("messages", [])

    await cl.Message(
        content=(
            "Welcome to **ChemPath** — your drug screening assistant.\n\n"
            "I can help you:\n"
            "- Screen compounds for protein targets (e.g., *Find EGFR inhibitors*)\n"
            "- Compute Pareto-optimal trade-offs (e.g., *What's the best balance for BRAF?*)\n"
            "- Run sensitivity analysis (e.g., *How robust are ALK rankings?*)\n"
            "- Look up compound details (e.g., *Tell me about Sirolimus*)\n\n"
            "What would you like to explore?"
        )
    ).send()


@cl.on_message
async def on_message(message: cl.Message):
    """Handle incoming user messages with Claude tool use."""
    client = cl.user_session.get("client")
    messages = cl.user_session.get("messages")

    messages.append({"role": "user", "content": message.content})

    msg = cl.Message(content="")
    await msg.send()

    # Agentic loop
    while True:
        response = await client.messages.create(
            model=MODEL,
            max_tokens=4096,
            system=SYSTEM_PROMPT,
            tools=TOOL_DEFINITIONS,
            messages=messages,
        )

        messages.append({"role": "assistant", "content": response.content})

        if response.stop_reason == "end_turn":
            text_parts = [
                block.text for block in response.content
                if hasattr(block, "text")
            ]
            msg.content = "\n".join(text_parts)
            await msg.update()
            break

        if response.stop_reason == "tool_use":
            tool_results = []
            for block in response.content:
                if block.type == "tool_use":
                    # Show tool call in UI
                    async with cl.Step(name=block.name, type="tool") as step:
                        step.input = block.input
                        result = executor.execute(block.name, block.input)
                        step.output = result
                        tool_results.append({
                            "type": "tool_result",
                            "tool_use_id": block.id,
                            "content": result,
                        })
            messages.append({"role": "user", "content": tool_results})
            continue

        # Unexpected stop
        msg.content = "(Unexpected response)"
        await msg.update()
        break

    cl.user_session.set("messages", messages)
