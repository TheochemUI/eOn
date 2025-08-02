# Tracking changes

`towncrier` is used to automate the generation and maintenance of the
`CHANGELOG.md` file. This ensures that all changes are properly attributed and
formatted, creating a consistent and readable history of the project's
evolution.

## Usage

The core workflow revolves around creating small "news fragments" for each
change. These fragments are then collected by towncrier to build the final
changelog entry for a release.

## Creating News Fragments

For each pull request that introduces a user-facing or noteworthy change, you
must add a news fragment.

- Navigate to the fragments directory: All fragments are stored in `docs/newsfragments/`.

- Create a new file: The filename must follow the format `PR_NUMBER.KEY.TYPE`.

- PR_NUMBER: The number of the corresponding pull request. If the change doesn't
  have a clear pull request, use a `+` symbol.

- TYPE: The category of the change. See the table below for available types.

- Write the fragment content: The file should contain a brief, past-tense
  description of the change, formatted in Markdown. The content will be rendered
  as a bullet point under the appropriate heading.

### Example

```
A fix for issue #456 is submitted in pull request #457. A new file would be
created at docs/newsfragments/457.fixed.md with the content:

Fixed a memory leak that occurred during long-running simulations.
```

## Fragment Types

The TYPE in the filename determines which section of the changelog the entry
will appear under. The following types are configured for this project:

| Type         | Heading    | Description                |
|--------------|------------|----------------------------|
| `added`      | Added      | New feature                |
| `changed`    | Changed    | Changes to features        |
| `deprecated` | Deprecated | Soon to be removed         |
| `removed`    | Removed    | Features which are gone    |
| `fixed`      | Fixed      | Bug fixes                  |
| `security`   | Security   | Vulnerability              |
| `dev`        | Developer  | Changes to the dev process |


## Generating the Changelog

When preparing for a new release, the towncrier command is used to assemble all
fragments into the `CHANGELOG.md` file.

Ensure your git working directory is clean before running the build command, as
`towncrier` will delete the fragment files it processes.

### Draft the Release

To preview the generated changelog notes without modifying any files, run the
build command with the `--draft` flag. This is useful for verifying that all
fragments are correctly formatted and included.

```
towncrier build --draft --version 0.1.0
```

### Build the Final Changelog

Once you are satisfied with the draft, run the command again without the `--draft`
flag. This will prepend the new release notes to `CHANGELOG.md` and delete the
consumed news fragment files.

```
towncrier build --version 0.1.0 --date "$(date -u +%Y-%m-%d)"
```
