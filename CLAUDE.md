# CLAUDE.md - AI Assistant Guide for ERS Repository

## 📋 Repository Overview

**Repository Name:** ERS
**Purpose:** [To be defined as project develops]
**Status:** Initial setup phase

This document serves as a comprehensive guide for AI assistants (like Claude) working on this codebase. It contains essential information about the project structure, conventions, and workflows to follow.

---

## 🏗️ Codebase Structure

```
ERS/
├── .git/                    # Git version control
├── CLAUDE.md               # This file - AI assistant guide
└── [Additional structure TBD]
```

### Directory Organization

**As this is a new repository, the structure will evolve. When adding directories, follow these conventions:**

- `/src` - Source code
- `/tests` or `/test` - Test files
- `/docs` - Documentation
- `/config` - Configuration files
- `/scripts` - Utility scripts
- `/public` or `/static` - Static assets (if web application)
- `/lib` or `/utils` - Shared utilities and libraries

---

## 🔧 Technology Stack

**To be determined based on initial commits. Update this section when technologies are chosen.**

Potential categories to document:
- **Language(s):**
- **Framework(s):**
- **Database:**
- **Build Tools:**
- **Testing:**
- **Deployment:**

---

## 🚀 Development Workflows

### Git Workflow

**Branch Naming Convention:**
- Feature branches: `claude/claude-md-[session-id]` (for Claude AI sessions)
- Feature branches: `feature/[feature-name]`
- Bug fixes: `fix/[bug-description]`
- Hotfixes: `hotfix/[issue]`

**Commit Message Format:**
```
[Type]: Brief description

Detailed explanation if needed
- Bullet points for multiple changes
- Keep it clear and concise
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `refactor`: Code refactoring
- `test`: Adding or updating tests
- `chore`: Maintenance tasks
- `style`: Code style changes (formatting, etc.)

### Pull Request Process

1. Create feature branch from main
2. Make changes with clear, atomic commits
3. Push to remote: `git push -u origin <branch-name>`
4. Create PR with descriptive title and summary
5. Wait for review and approval
6. Merge when approved

---

## 🤖 Key Conventions for AI Assistants

### Code Quality Standards

1. **Write Clean, Readable Code**
   - Use meaningful variable and function names
   - Keep functions small and focused (single responsibility)
   - Comment complex logic, but prefer self-documenting code
   - Follow DRY (Don't Repeat Yourself) principle

2. **Error Handling**
   - Always handle errors appropriately
   - Use try-catch blocks where necessary
   - Provide meaningful error messages
   - Log errors for debugging

3. **Testing**
   - Write tests for new features
   - Ensure existing tests pass before committing
   - Aim for meaningful test coverage (quality over quantity)

4. **Security**
   - Never commit secrets, API keys, or credentials
   - Use environment variables for sensitive data
   - Validate and sanitize all inputs
   - Follow security best practices (OWASP guidelines)

### File Operations

**Prefer Editing Over Creating:**
- Always check if a file exists before creating a new one
- Use Edit tool for modifications to existing files
- Only create new files when absolutely necessary

**Reading Before Writing:**
- Always read a file before modifying it
- Understand the existing code structure
- Maintain consistency with existing patterns

### Communication

**Be Concise:**
- Keep explanations clear and brief
- Focus on what changed and why
- Use bullet points for clarity

**Code References:**
- Reference specific files and line numbers
- Format: `file_path:line_number`
- Example: `src/services/auth.js:45`

---

## 📝 Common Development Tasks

### Starting a New Feature

```bash
# Create and switch to feature branch
git checkout -b feature/your-feature-name

# Make your changes
# ... code, test, commit ...

# Push to remote
git push -u origin feature/your-feature-name
```

### Running Tests

```bash
# Update this section when test framework is established
# Example: npm test, pytest, go test, etc.
```

### Building the Project

```bash
# Update this section when build process is established
# Example: npm run build, make, go build, etc.
```

### Running Locally

```bash
# Update this section when development server is established
# Example: npm start, python app.py, etc.
```

---

## 🎯 Project-Specific Guidelines

### Naming Conventions

**Update as patterns emerge:**
- **Variables:** camelCase or snake_case (TBD)
- **Functions:** camelCase or snake_case (TBD)
- **Classes:** PascalCase
- **Constants:** UPPER_SNAKE_CASE
- **Files:** kebab-case or snake_case (TBD)

### Code Style

**Will be defined based on chosen language/framework:**
- Indentation: [spaces/tabs] [number]
- Line length: [max characters]
- Quotes: [single/double]
- Semicolons: [required/optional]

Consider using:
- Linter configuration (ESLint, Pylint, etc.)
- Formatter (Prettier, Black, gofmt, etc.)
- EditorConfig for consistency

---

## 🔍 Important Notes for AI Assistants

### Before Making Changes

1. **Read the relevant files first** - Never propose changes to code you haven't read
2. **Understand the context** - Know how your changes fit into the larger system
3. **Check for existing patterns** - Follow established conventions in the codebase
4. **Search for related code** - Ensure consistency across similar functionality

### When Implementing Features

1. **Avoid over-engineering** - Keep solutions simple and focused
2. **Don't add unrequested features** - Stick to the requirements
3. **Minimal necessary changes** - Don't refactor unrelated code
4. **Security first** - Watch for vulnerabilities (XSS, SQL injection, etc.)

### Git Operations Best Practices

1. **Branch naming:** Must start with `claude/` for AI sessions
2. **Always use:** `git push -u origin <branch-name>`
3. **Retry on network errors:** Up to 4 times with exponential backoff (2s, 4s, 8s, 16s)
4. **Never skip hooks:** Unless explicitly requested
5. **Verify before force operations:** Never force push to main/master

### Commit Strategy

1. **Clear, descriptive messages** - Explain the "why" not just the "what"
2. **Atomic commits** - One logical change per commit
3. **Follow the commit format** - Use conventional commit types
4. **Review before committing** - Check `git status` and `git diff`

---

## 📚 Resources and Documentation

### External Documentation Links
[Add relevant links as project develops]

- Project Management: [Link TBD]
- API Documentation: [Link TBD]
- Design System: [Link TBD]
- Deployment Docs: [Link TBD]

### Internal Documentation
[Add as documentation is created]

- Architecture Decision Records (ADRs): [Location TBD]
- API Specifications: [Location TBD]
- Database Schema: [Location TBD]

---

## 🔄 Updating This Document

This document should be updated whenever:
- New technologies or frameworks are added
- Development workflows change
- New conventions are established
- Project structure evolves significantly

**Last Updated:** 2025-11-28
**Version:** 1.0.0 (Initial)

---

## 🎓 Quick Reference for Common Scenarios

### Scenario 1: User asks to add a new feature
1. Read relevant existing code first
2. Plan the implementation (use TodoWrite if complex)
3. Make minimal necessary changes
4. Test the changes
5. Commit with clear message
6. Push to feature branch

### Scenario 2: User reports a bug
1. Locate the problematic code
2. Understand the root cause
3. Fix the issue without over-refactoring
4. Add tests if appropriate
5. Commit and push

### Scenario 3: User asks about codebase
1. Use Explore agent for broad questions
2. Use Grep/Glob for specific searches
3. Read relevant files
4. Provide clear, concise explanation with file references

### Scenario 4: Creating a Pull Request
1. Ensure all changes are committed
2. Push to remote branch
3. Review git diff to understand full scope
4. Create PR with summary of all commits (not just latest)
5. Include test plan in PR description

---

**Remember:** This is a living document. As the ERS project grows and evolves, keep this guide updated to reflect the current state of the codebase and best practices.
