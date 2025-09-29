# Git Version Control Guidelines 🔀

Welcome to the PANORAMA git workflow guide! Whether you're a core maintainer or a first-time contributor, this document
will help you navigate our development process. We're excited to have you here and want to make contributing as smooth
as possible.

```{note}
**New to contributing?** Don't worry! We're here to help. If anything in this guide is unclear or you get stuck, just reach out. We've all been there, and we're happy to guide you through the process.
```

## Understanding Our Branch Structure 🌿

PANORAMA uses a two-branch model to separate stable releases from active development:

- **`main`** - The stable, production-ready version. This is what users download and run.
- **`dev`** - The latest development version where all new work comes together. This is your starting point!

Think of `dev` as the collaborative workspace where features come together before they're polished and released to
`main`.

### Branch Naming Conventions

When creating your feature branch from `dev`, use descriptive names with these prefixes:

- `feature/` - New functionality or enhancements
- `bugfix/` - Bug fixes
- `docs/` - Documentation updates
- `refactor/` - Code improvements without changing functionality
- `test/` - Adding or updating tests
- `hotfix/` - Urgent fixes (usually branched from `main`)

**Examples:**

```bash
feature/add-ai-algorithms
bugfix/fix-annotation-memory-leak
docs/improve-installation-guide
refactor/optimize-clustering-algorithm
```

### Typical Branch Lifecycle

Here's how a feature makes its way into PANORAMA:

```bash
# 1. Start from the latest dev
git checkout dev
git pull origin dev

# 2. Create your feature branch
git checkout -b feature/my-awesome-feature

# 3. Develop, test, and document your changes
# ... write code, add tests, update docs ...

# 4. Push your branch
git push -u origin feature/my-awesome-feature

# 5. Create a Pull Request to merge into dev
# Your PR gets reviewed and merged into dev

# 6. Eventually, accumulated changes in dev get released
# A maintainer merges dev → main with a new minor version
```

After your PR is merged, feel free to delete your feature branch - it's served its purpose!

```bash
git branch -d feature/my-awesome-feature
git push origin --delete feature/my-awesome-feature
```

```{important}
The reviewer can also delete your branch after merging direclty from GitHub.
```

## Versioning and Releases 🏷️

PANORAMA follows semantic versioning: `Major.Minor.Patch`

- **Patch** (e.g., 1.2.**3**) - Each Pull Request increments the patch number. Bug fixes, small improvements,
  documentation updates.
- **Minor** (e.g., 1.**3**.0) - When `dev` is merged into `main`, the minor version bumps. New features, significant
  improvements.
- **Major** (e.g., **2**.0.0) - Breaking changes or major architectural updates. Rare but important!

**What this means for your PR:**
When you create a pull request, you'll need to bump the patch version in the `VERSION` file. Don't worry - reviewers
will help you if you forget!

```bash
# Current version: 1.2.5
# Your PR changes it to: 1.2.6
echo "1.2.6" > VERSION
git add VERSION
git commit -m "Bump version to 1.2.6"
```

## Writing Good Commits 📝

Good commit messages are like good lab notes - they help everyone (including future you) understand what happened and
why.

### The Basic Format

**One-line commits** for simple changes:

```bash
git commit -m "Fix off-by-one error in gene counting"
git commit -m "Add validation for empty gene families"
git commit -m "Update installation instructions for conda"
```

**Multi-line commits** when you need to explain more:

```bash
git commit -m "Optimize system clustering for large datasets

The previous implementation used nested loops that didn't scale well.
This commit introduces vectorized operations and caching that reduce
runtime from 2 hours to 25 minutes on 10k+ genome datasets.

Tested on: E. coli dataset (15k genomes), Staphylococcus (8k genomes)"
```

### Commit Message Tips

- **Use imperative mood**: "Add feature" not "Added feature"
- **Be specific**: "Fix memory leak in clustering" beats "Fix bug"
- **Keep first line under 50 characters** when possible
- **Explain the 'why' not the 'how'** - code shows how, commits explain why
- **Make atomic commits** - one logical change per commit

```{tip}
 If you find yourself using "and" in a commit message, you might want to split it into multiple commits!
```

### Small, Focused Commits

Break your work into digestible pieces:

```bash
# Good: Three clear, focused commits
git commit -m "Add merge method to System class"
git commit -m "Add unit tests for System.merge()"  
git commit -m "Document System.merge() in API reference"

# Less ideal: One big blob
git commit -m "Add merge feature with tests and docs"
```

This makes it easier to review, debug, and potentially revert changes if needed.

## Before You Push: The Checklist ✅

We've all pushed code and then immediately realized something was wrong. This checklist helps catch issues before they
become embarrassing! 😅

### 1. Run the Tests

```bash
# Quick check
pytest

# Full check with coverage (recommended)
pytest --cov=panorama

# Just test what you changed
pytest tests/test_my_feature.py
```

All tests should pass. If something fails, fix it before pushing. Your future self will thank you!

### 2. Format with Black

We use Black to keep code style consistent. No more debates about spaces and brackets!

```bash
# Format everything
black panorama/ tests/

# Check what would change (without modifying files)
black --check panorama/ tests/
```

Black makes code reviews smoother since we're focusing on logic, not style.

### 3. Linting with Flake8
[flake8](https://flake8.pycqa.org/) catches potential bugs and style issues:

```shell
# Check the entire project
flake8 panorama/ tests/

# Check specific files
flake8 panorama/systems/system.py
```

Fix the issues flake8 reports before pushing. Most are quick fixes!

### 4. Update Documentation

Documentation is code too! If you:

- **Added a new feature** → Update user documentation
- **Changed public APIs** → Update API reference
- **Added/modified functions** → Write/update docstrings

More information on how to write good documentation can be found in the 
["how to build the documentation"](buildDoc.md).

### 5. Review Your Own Changes

Before asking others to review, review yourself:

```bash
# What changed compared to dev?
git diff origin/dev

# Check your commit history
git log origin/dev..HEAD --oneline

# Make sure you didn't leave any debug code
grep -r "print(" panorama/  # Just an example!
```

### 6. Update the VERSION File

Don't forget to bump the patch version! See the [Versioning section](#versioning-and-releases-🏷️) above.

## Creating a Pull Request 🔄

You've done the hard work - now let's get it merged!

### PR Title and Description

Write a clear title that summarizes your change:

```
✅ Good titles:
- Add support for detecting CRISPR-Cas systems
- Fix memory leak in pangenome clustering  
- Update testing guidelines with pytest examples

❌ Less helpful titles:
- Update
- Fix bug
- Changes
```

Use the PR description to provide context:

```markdown
## What does this PR do?

Adds a new module for detecting CRISPR-Cas systems in bacterial genomes using HMM profiles.

## Why is this needed?

Many users requested CRISPR detection capabilities. This addresses issue #123.

## What changed?

- Added `panorama/crispr/detector.py` with detection logic
- Added HMM profiles for Cas proteins
- Added unit tests (coverage: 95%)
- Updated user documentation with examples

## How to test?

pytest tests/test_crispr.py
panorama detect_crispr --pangenome examples/ecoli.h5

## Related Issues

Closes #123
Related to #456

## Version

Updated VERSION from 1.2.5 → 1.2.6
```

### Linking Issues

If your PR addresses an issue, link it! Use these keywords in your PR description:
- `Fixes #123` - Automatically closes issue when PR merges
- `Closes #123` - Same as above
- `Related to #123` - Links without closing

## Code Review: A Collaborative Process 👥

Code review isn't about finding faults - it's about improving code together and sharing knowledge!

### For Authors (You!)

- **Be responsive** to feedback, but don't feel pressured to accept every suggestion
- **Ask questions** if feedback is unclear - reviewers aren't always right!
- **Don't take it personally** - comments are about code, never about you
- **Mark conversations as resolved** once you've addressed them

```{hint}
**Remember:** Even experienced developers get feedback on their PRs. It's part of the process and makes everyone better!
```

### For Reviewers

All PRs need at least one maintainer review before merging. If you're reviewing:

**Check these things:**
1. Does the code work as intended?
2. Are there tests? Do they cover edge cases?
3. Is documentation updated?
4. Is the code readable and maintainable?
5. Are there any performance or security concerns?
6. Is the VERSION file updated?

**How to give good feedback:**
```markdown
✅ Constructive:
"This could cause issues with empty datasets. Could we add a check for len(dataset) > 0?"

"Nice solution! FYI, there's a built-in function that does this: `collections.Counter`"

"I'm curious why we're using approach X instead of Y here - could you explain the reasoning?"

❌ Less helpful:
"This is wrong."
"Bad code."
"Why would you do it this way?"
```

Be kind, be specific, be constructive. We're all learning together!

## Keeping Your Branch Updated 🔄

While you're working on your feature, `dev` keeps moving forward. Here's how to stay in sync:

```bash
# Update your local dev
git checkout dev
git pull origin dev

# Bring those changes into your feature branch
git checkout feature/my-feature
git rebase dev

# If you've already pushed your branch
git push --force-with-lease
```

**Rebase vs Merge:** We generally prefer rebasing for a cleaner history, but merging is fine too. Use what you're
comfortable with!

### Handling Merge Conflicts

Conflicts happen - they're not a failure, just Git asking for your help to combine changes.

```bash
# Start the rebase
git rebase dev
# Git pauses on conflicts

# Open conflicting files and look for markers:
<<<<<<< HEAD
Your code
=======
Their code
>>>>>>> dev

# Edit to keep what you want, then:
git add resolved_file.py
git rebase --continue

# If things get messy, you can always abort and try again
git rebase --abort
```

**Stuck on conflicts?** Don't hesitate to ask for help! Ping a maintainer or open a draft PR and explain where you're
stuck.

## Common Workflows 🔄

### Adding a New Feature

```bash
# 1. Start fresh from dev
git checkout dev
git pull origin dev
git checkout -b feature/awesome-new-thing

# 2. Build your feature
# ... code code code ...
git add panorama/awesome/feature.py
git commit -m "Add awesome new feature"

# ... test test test ...
git add tests/test_awesome.py
git commit -m "Add tests for awesome feature"

# ... document document document ...
git add docs/user/awesome_guide.md
git commit -m "Document awesome feature usage"

# 3. Pre-push checks
pytest
black panorama/ tests/

# 4. Update version
echo "1.2.6" > VERSION  # Assuming current is 1.2.5
git add VERSION
git commit -m "Bump version to 1.2.6"

# 5. Push and create PR
git push -u origin feature/awesome-new-thing
# Go to GitHub and create your PR!
```

### Fixing a Bug

```bash
# 1. Create bugfix branch from dev
git checkout dev
git pull origin dev
git checkout -b bugfix/fix-that-annoying-thing

# 2. Fix it
git add panorama/problematic_module.py
git commit -m "Fix off-by-one error in gene counting

The previous implementation incorrectly counted the last gene in each
family. This fix ensures all genes are counted correctly.

Fixes #234"

# 3. Add a test so it doesn't come back
git add tests/test_gene_counting.py
git commit -m "Add regression test for gene counting bug"

# 4. Version bump and push
echo "1.2.7" > VERSION
git add VERSION
git commit -m "Bump version to 1.2.7"

pytest
black panorama/ tests/
git push -u origin bugfix/fix-that-annoying-thing
```

### Updating Documentation Only

```bash
git checkout -b docs/improve-installation-guide

# Make your changes
git add docs/user/installation.md
git commit -m "Clarify conda installation steps for macOS"

# Documentation-only changes still need version bumps
echo "1.2.8" > VERSION
git add VERSION  
git commit -m "Bump version to 1.2.8"

git push -u origin docs/improve-installation-guide
```

## Useful Git Commands 🛠️

Some handy commands to make your life easier:

```bash
# Beautiful commit history
git log --oneline --graph --all

# What changed in a specific commit?
git show abc123

# Temporarily save changes without committing
git stash
git stash pop  # Get them back

# Oops, need to change the last commit message?
git commit --amend -m "Better message"

# Interactive rebase to clean up commits before pushing
git rebase -i HEAD~3  # Last 3 commits

# Find which commit introduced a bug
git bisect start
```

## Getting Help 🆘

**Stuck? Confused? Not sure what to do?**

We're here to help! Seriously. Here are your options:

1. **Open a draft PR** - Get early feedback even if it's not ready
2. **Ask in discussions** - No question is too basic
3. **Tag a maintainer** - We're friendly, promise!
4. **Check the documentation** - But don't spend hours lost in docs
5. **Open an issue** - "Help wanted" tags are specifically for questions

```{tip}
It's much better to ask a "dumb" question (there are no dumb questions!) than to struggle in silence or
make a mistake that's hard to fix later.
```

## For External Contributors 🌍

We genuinely love external contributions! Here's what you should know:

**First time contributing?**

- Start with issues labeled `good first issue` or `help wanted`
- Don't be intimidated by the process - we'll guide you through
- It's okay if your first PR isn't perfect - that's what reviews are for!

**Not sure if your idea fits?**

- Open an issue first to discuss it
- We're very open to new ideas and improvements
- Even if we don't merge something, the discussion is valuable

**Need help with the workflow?**

- Ask! We're happy to explain git, pytest, documentation, whatever
- We can help you get your dev environment set up
- Video calls are fine if you're really stuck

Remember: Every expert was once a beginner. We all started somewhere, and we want to support your journey! 🚀

## Quick Reference Card 📋

```bash
# Start new feature
git checkout dev && git pull origin dev
git checkout -b feature/my-feature

# Before pushing
pytest
black panorama/ tests/
# Update VERSION file
git push -u origin feature/my-feature

# Keep branch updated  
git checkout dev && git pull origin dev
git checkout feature/my-feature
git rebase dev

# After PR is merged
git branch -d feature/my-feature
git push origin --delete feature/my-feature
```

Happy coding! 🎉