"""Backticks inside a ``$( ... )`` command substitution in the gate's own shell scripts.

**The defect.** Inside a double-quoted command substitution, a backtick opens a *nested* substitution.
So a markdown-style code span written in a **comment** — natural in this repository, where the shell
scripts carry long prose — makes the script emit ``syntax error: unexpected end of file`` at runtime.

**Why a check and not care.** ``bash -n`` does not catch it. Only running the script does, and the
script still produces most of its output, so the error scrolls past. It reached committed ``main`` on
14.08.2026 (``60ba06c``) and survived until someone read the gate's output line by line.

**Why this is a separate file.** It was first written inline in ``verify_main.sh`` as a heredoc inside
``B=$( ... )`` — and the Python needed a backtick to search for one, which broke the enclosing
substitution exactly as described above. The check reproduced its own defect twice before landing
here. A helper file has no enclosing substitution and cannot.

**Backticks in ordinary comments outside a substitution are harmless and are not flagged.**

**Verified by fault injection**, which is this gate's standard: a backtick was inserted inside the
multi-line ``M=$(.venv/bin/python -c " ... ")`` block and the check must report it. The first version
did **not** — its regex required the closing parenthesis to start a line, so it skipped precisely the
multi-line substitutions the defect lives in, and reported 0. That is why this scans with matched
parentheses instead.

⚠️ **Known limits.** It is a bracket scanner, not a shell parser: a literal ``)`` inside a quoted
string within a substitution will end the span early, so a backtick after that point in the same
substitution is missed. It reports how many spans it examined, so a collapse to zero is visible --
the failure mode that had the first version reporting 0 against a true 11 backticks.

Usage:  shell_safety.py <gate_dir>
"""
import sys
from pathlib import Path

BACKTICK = chr(96)


def substitutions(text: str):
    """Yield ``(start_index, body)`` for every ``$( ... )`` span, matching parentheses properly.

    A regex cannot do this: the first version of this file used one that required the closing
    parenthesis to start a line, and it therefore skipped exactly the multi-line ``X=$(python -c "
    ... ")`` blocks where the defect lives. Fault injection caught that -- the checker reported zero
    against an injected backtick -- which is why this is a scanner and why the injection test is
    described in the module docstring.
    """
    i, n = 0, len(text)
    while i < n - 1:
        if text[i] == "$" and text[i + 1] == "(":
            depth, j = 1, i + 2
            start = j
            while j < n and depth:
                c = text[j]
                if c == "(":
                    depth += 1
                elif c == ")":
                    depth -= 1
                j += 1
            yield i, text[start:j - 1]
            i = j
        else:
            i += 1


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    files = sorted(root.glob("*.sh"))
    spans = 0
    findings = []
    for f in files:
        text = f.read_text(errors="ignore")
        for pos, body in substitutions(text):
            spans += 1
            if BACKTICK in body:
                line = text[:pos].count("\n") + 1
                findings.append(f"   BACKTICK IN SUBSTITUTION  {f.name}:{line}")
    print(f"  gate shell: {len(files)} script(s), {spans} command substitution(s), "
          f"{len(findings)} containing a backtick")
    for x in findings:
        print(x)
    if not files or spans == 0:
        print("   ^ BLOCKER: nothing examined -- '0 found' would be indistinguishable from "
              "'0 looked at'")
        return 2
    return 1 if findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
