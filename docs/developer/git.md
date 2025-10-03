(git-guide)=
# Git Version Control Guidelines 🗃️

Welcome to the PANORAMA git workflow guide! Whether you're a core maintainer or a first-time contributor, this document
will help you navigate our development process. We're excited to have you here and want to make contributing as smooth
as possible.

## Understanding Our Branch Structure 🌿

PANORAMA uses a two-branch model to separate stable releases from active development:

- **`main`** - The stable, production-ready version. This is what users download and run.
- **`dev`** - The latest development version where all new work comes together. This is your starting point!

Think of `dev` as the collaborative workspace where features come together before they're polished and released to
`main`.

### Branch Naming Conventions

When creating your feature branch from `dev`, use descriptive names.

**Examples:**

```bash
add-ai-algorithms
fix-annotation-memory-leak
improve-installation-guide
optimize-clustering-algorithm
```

```{tip}
You can add a prefixe to your branch name to indicate its purpose:

- `feature/` - New functionality or enhancements
- `bugfix/` - Bug fixes
- `docs/` - Documentation updates
- `refactor/` - Code improvements without changing functionality
- `test/` - Adding or updating tests
- `hotfix/` - Urgent fixes (usually branched from `main`)
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

VERSION will be updated by the maintainers when merging PRs.
But you can also update it manually to make his life easier (he will appreciate it):

```bash
# Current version: 1.2.5
# Your PR changes it to: 1.2.6
echo "1.2.6" > VERSION
git add VERSION
git commit -m "Bump version to 1.2.6"
```


## Creating a Pull Request 🔄

You've done the hard work - now let's get it merged!

### PR Title and Description

Write a clear title that summarizes your change and use the PR template to provide context:

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

```{important}
If your PR addresses an issue, link it! Use these keywords in your PR description:
- `Fixes #123` - Automatically closes issue when PR merges
- `Closes #123` - Same as above
- `Related to #123` - Links without closing
```

## Keeping Your Branch Updated

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

## For Contributors 🌍

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
