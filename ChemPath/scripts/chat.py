"""
ChemPath CLI Chat — interactive terminal interface.
Requires ANTHROPIC_API_KEY environment variable.
"""

import sys
from chempath.agent.engine import ChemPathAgent


def main():
    print("=" * 70)
    print(" ChemPath — Drug Screening Chat Assistant")
    print("=" * 70)
    print("  Type your question about drug compounds and targets.")
    print("  Commands: /reset (clear history), /quit (exit)")
    print("=" * 70)

    try:
        agent = ChemPathAgent()
    except ImportError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        if "API key" in str(e) or "api_key" in str(e).lower():
            print("Error: Set ANTHROPIC_API_KEY environment variable.")
            print("  export ANTHROPIC_API_KEY=sk-ant-...")
            sys.exit(1)
        raise

    while True:
        try:
            user_input = input("\nYou: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGoodbye!")
            break

        if not user_input:
            continue
        if user_input == "/quit":
            print("Goodbye!")
            break
        if user_input == "/reset":
            agent.reset()
            print("  [Conversation reset]")
            continue

        print("\nChemPath: ", end="", flush=True)
        try:
            for chunk in agent.chat_stream(user_input):
                print(chunk, end="", flush=True)
            print()
        except Exception as e:
            print(f"\n  [Error: {e}]")


if __name__ == "__main__":
    main()
