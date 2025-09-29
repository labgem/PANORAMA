# Reporting Issues and Suggesting Features 💬

Your feedback helps make PANORAMA better! Whether you've found a bug, have an idea for a new feature, or want to suggest an improvement, we want to hear from you. This guide will show you how to report issues and make suggestions effectively.

## Reporting a Bug 🐛

Found something that doesn't work as expected? Here's how to report it:

### Before Reporting

1. **Check if it's already reported** - Search the [issue tracker](https://github.com/labgem/PANORAMA/issues) to see if someone else has reported the same issue
2. **Make sure it's a bug** - Is this a problem with PANORAMA, or an issue with your input data or environment?
3. **Try to reproduce it** - Can you make the bug happen again with the same steps?

### Creating a Bug Report

Go to the [PANORAMA issue tracker](https://github.com/labgem/PANORAMA/issues/new) and use this template:

---
```markdown
## Bug Description
A clear description of what the bug is.

## Steps to Reproduce
1. Run this command: `panorama ...`
2. With these input files: `...`
3. See error

## Expected Behavior
What you expected to happen.

## Actual Behavior
What actually happened. Include complete error messages.

## Environment
- PANORAMA version: [run `panorama --version`]
- Python version: [run `python --version`]
- Operating System: [e.g., Ubuntu 22.04, macOS 13, Windows 11]
- Installation method: [pip, conda, from source]

## Input Data (if applicable)
- Number of genomes:
- Pangenome size (number of gene families):
- Any specific characteristics of your data:

## Error Message
Paste the complete error message here
Include the full traceback if available

## Additional Context
Any other information that might help us understand or reproduce the issue.

## Possible Solution (optional)
If you have ideas about what might be causing this or how to fix it.
```
---

## Suggesting a New Feature 💡

Have an idea to improve PANORAMA? We'd love to hear it!

### Before Suggesting

1. **Check existing features** - Make sure the feature doesn't already exist
2. **Search the issues** - Someone might have already suggested something similar
3. **Think about the use case** - How would this feature be useful to you and others?

### Creating a Feature Request

Go to the [PANORAMA issue tracker](https://github.com/labgem/PANORAMA/issues/new) and use this template:

---
```markdown
## Feature Description
A clear description of the feature you'd like to see.

## Use Case
Describe the problem this feature would solve or the workflow it would improve.

**Example scenario:**
"When I'm analyzing large pangenomes with >10k genomes, I need to..."

## Proposed Solution
How do you envision this feature working?

**Example:**
`panorama new_command --option value`

## Alternatives Considered
Have you thought of any alternative solutions or workarounds?

## Additional Context
Any other information, mockups, examples from other tools, papers, etc.

## Would you be willing to contribute?
- [ ] Yes, I'd like to work on this
- [ ] Yes, with guidance
- [ ] No, but I'd be happy to test it
- [ ] Just suggesting the idea
```
---

## Enhancement Suggestion ✨

Some idea to enhance PANORAMA? Tell us about it!
```markdown
## What to Improve
Which part of PANORAMA could be better?

**Area:** [Documentation / Performance / Usability / Output / etc.]

## Current Situation
How does it work now, and what's the issue?

## Suggested Improvement
What would make it better?

## Expected Benefits
Who would benefit and how?

## Examples (optional)
Examples from other tools or how you imagine it could work.
```
---

## Asking Questions or Sharing Ideas 💭

Not sure if something is a bug or feature request? Just want to discuss an idea? Use GitHub Discussions!

Go to [GitHub Discussions](https://github.com/labgem/PANORAMA/discussions) and share:

- Questions about how PANORAMA works
- Ideas you want to explore before making a formal feature request
- General feedback or suggestions
- Use cases and workflows you'd like to share

---
```markdown
## Topic
A brief, clear title for your discussion.

## Context
What prompted this question or idea?

## Question/Idea
Your question or idea in detail.

**Example questions:**
- How does PANORAMA handle...?
- What's the best way to...?
- Has anyone tried...?

**Example ideas:**
- What if we could...?
- I've been thinking about...?
- Would it make sense to...?

## Your Thoughts
Any initial thoughts, attempts, or research you've done.

## Relevant Information (optional)
- Related papers or tools
- Example data or use cases
- Potential challenges or considerations
```
---

## What Happens Next? 🔄

After you submit an issue or discussion:

1. **We'll respond** - Usually within a few days
2. **We might ask questions** - To better understand your issue or idea
3. **We'll label it** - To help organize and prioritize
4. **We'll discuss solutions** - For bugs, or feasibility for features
5. **Someone might work on it** - Could be a maintainer or a contributor (maybe you!)

```{note}
Be patient - we're a small team, but we read every issue and discussion!
```

## Tips for Good Reports 📋

### Do:
- Be specific and detailed
- Include examples and use cases
- Provide complete error messages
- Be respectful and constructive
- Follow up if we ask questions

### Don't:
- Report multiple unrelated issues in one report
- Demand immediate fixes or features
- Include sensitive or private data
- Duplicate existing issues

---

## Example Issues 📚

### Good Bug Report Example

```markdown
## Bug Description
The `panorama systems` command crashes with a KeyError when processing 
pangenomes that have gene families without functional annotations.

## Steps to Reproduce
1. Create a pangenome without running annotation first
2. Run: `panorama systems --pangenomes test.h5 --models defense_models.yml`
3. Error occurs during system detection

## Expected Behavior
Should either skip unannotated families or show a helpful error message 
suggesting to run annotation first.

## Actual Behavior
    ```shell
    KeyError: 'annotation'
    Traceback (most recent call last):
      File "panorama/systems/detection.py", line 145, in detect_systems
        annotation = family.metadata['annotation']
    KeyError: 'annotation'
    ```

## Environment
- PANORAMA version: 1.2.3
- Python version: 3.10.8
- Operating System: Ubuntu 22.04
- Installation: pip

## Input Data
- 500 genomes
- 12,000 gene families
- No annotation step performed

## Possible Solution
Check if annotation metadata exists before accessing it, or add a validation 
step that ensures annotation has been run before system detection.
```

### Good Feature Request Example

```markdown
## Feature Description
Add support for exporting system detection results in JSON format.

## Use Case
I'm integrating PANORAMA into an automated pipeline that processes results 
with Python scripts. Currently, I have to parse TSV files which is error-prone.
JSON would allow direct loading into Python dictionaries and easier integration 
with web APIs.

## Proposed Solution
Add a `--format` option to `panorama write_systems`:

    ```bash
    panorama write_systems --pangenomes data.h5 --output results/ --format json
    ```

Output structure:
    ```json
    {
      "systems": [
        {
          "id": "system_001",
          "model": "CRISPR-Cas",
          "organisms": ["org1", "org2"],
          "completeness": 0.95,
          "gene_families": [...]
        }
      ]
    }
    ```

## Alternatives Considered
- Parsing TSV files (current workaround, but fragile)
- Converting TSV to JSON with external tools (extra step)

## Additional Context
Many bioinformatics tools now offer JSON output (e.g., Prokka, eggNOG-mapper).
This would align PANORAMA with modern pipeline practices.

## Would you be willing to contribute?
- [x] Yes, with guidance - I could implement this with some pointers on where 
  to add the JSON export logic.
```

---

## Getting Help 🆘

If you're not sure how to report something or have questions:

- Check the [documentation](link-to-docs)
- Ask in [GitHub Discussions](https://github.com/labgem/PANORAMA/discussions)
- Open an issue anyway - we'll help you refine it!

Thank you for helping improve PANORAMA! Your feedback makes the tool better for everyone. 🙏
