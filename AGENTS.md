# AGENTS.md - Codex Full-Stack Dev Environment

This project includes pre-installed Codex skills (`.codex/skills/`) covering the full development lifecycle.

## Available Skills

| Skill | Trigger Scenarios | Path |
|-------|-------------------|------|
| `api-first-modular` | Frontend/backend development, cross-layer task decomposition, API design | `.codex/skills/api-first-modular/` |
| `code-debugger` | Bug fixing, performance tuning, incremental development | `.codex/skills/code-debugger/` |
| `debug-ui` | Frontend UI debugging - styling, interaction, rendering issues | `.codex/skills/debug-ui/` |
| `ai-spec` | Natural language requirements to precise technical spec translation | `.codex/skills/ai-spec/` |
| `ralph` | Autonomous development loop driven by PRD | `.codex/skills/ralph/` |

## Core Rules

- **API-First**: Every backend feature must be encapsulated as an independent API package (Implement -> Checkfix -> Encapsulate -> Expose API -> Document API). Frontend only consumes APIs per documentation - no business logic in the frontend.
- **Layer-scoped debugging**: Identify the bug's owning layer before making changes. Never apply cross-layer workarounds.
- **Ordered cross-layer execution**: Backend first -> API docs -> frontend consumption -> integration verification.

## Skills

A skill is a set of local instructions stored in a `SKILL.md` file. The list below includes name, description, and path so you can open each source when needed.

### Available skills

- ai-spec: Converts natural language requirements into production-grade technical specs and executable AI instructions. Use when requirements are unclear, architecture choices are needed, or a complete task list is required. (file: `E:/Development/Squidiff/.codex/skills/ai-spec/SKILL.md`)
- api-first-modular: API-first modular workflow. Backend features are packaged as independent APIs, frontend only calls APIs, and cross-layer work is split by API boundaries. Use for backend feature work, frontend API consumption, cross-layer decomposition, and bug fixing after layer ownership is identified. (file: `E:/Development/Squidiff/.codex/skills/api-first-modular/SKILL.md`)
- code-debugger: Deep-context debugging and incremental development workflow. Use for bug fixing, performance issues, variable tracing, and iterative extension of existing modules while maintaining `.debug/` records. (file: `E:/Development/Squidiff/.codex/skills/code-debugger/SKILL.md`)
- debug-ui: UI design and implementation workflow driven by product vibe. Converts high-level visual intent into precise CSS/design tokens and records design decisions (ADR) in `.debug/`. (file: `E:/Development/Squidiff/.codex/skills/debug-ui/SKILL.md`)
- ralph: PRD-driven autonomous agent loop. Converts Markdown PRD into `prd.json`, then iterates user stories with fresh agent runs until completion. (file: `E:/Development/Squidiff/.codex/skills/ralph/SKILL.md`)
- skill-creator: Guide for creating or updating skills that extend Codex capabilities with specialized workflows or integrations. (file: `C:/Users/DamnCheater/.codex/skills/.system/skill-creator/SKILL.md`)
- skill-installer: Installs skills into `$CODEX_HOME/skills` from curated lists or GitHub paths, including private repos. (file: `C:/Users/DamnCheater/.codex/skills/.system/skill-installer/SKILL.md`)

### How to use skills

- Discovery: The list above defines the skills available in this session. Skill bodies live at the listed paths.
- Trigger rules: If the user names a skill (with `$SkillName` or plain text), or if the request clearly matches a listed skill description, use that skill for the current turn. If multiple skills are named, use all relevant ones. Do not carry skills across turns unless re-mentioned.
- Missing/blocked: If a named skill is unavailable or unreadable, state it briefly and continue with the best fallback.
- How to use a skill (progressive disclosure):
  1) After selecting a skill, open its `SKILL.md` and read only what is needed.
  2) Resolve relative paths in `SKILL.md` against that skill directory first.
  3) If extra folders (for example, `references/`) are mentioned, load only specific files needed for the task.
  4) If `scripts/` exist, prefer running or patching them over retyping large code blocks.
  5) Reuse assets/templates when available.
- Coordination and sequencing:
  - If multiple skills apply, choose the minimal set and state the execution order.
  - Announce which skills are used and why in one short line. If skipping an obvious skill, say why.
- Context hygiene:
  - Keep context small. Summarize long sections and load extra files only when required.
  - Avoid deep reference chasing. Open files directly linked from `SKILL.md` unless blocked.
  - When variants exist (framework/provider/domain), choose only the relevant references and note the choice.
- Safety and fallback: If a skill cannot be cleanly applied (missing files, unclear instructions), state the issue and continue with the next-best approach.
