# Migrating from SVN

```{note}
This section should only be relevant until the `2.x` release at which point the SVN repository shall be completely removed.
```

The basic idea is to use `git svn`[^1], following the standard [git book
approach](https://git-scm.com/book/en/v2/Git-and-Other-Systems-Migrating-to-Git).

## Getting the SVN sources

For the first part we don't need **every** commit, so a checkout will work over
`clone`.

```{code-block} bash
svn co https://theory.cm.utexas.edu/svn/eon eonSVNonly
```

We will use this to determine author information.

```{code-block} bash
cd eonSVNonly
svn log > loggy
svn log --xml --quiet | grep author | sort -u | \
  perl -pe 's/.*>(.*?)<.*/$1 = /' > ../authors.txt
```

Now, the `authors.txt` file can be edited to map to existing GitHub users. It is
also useful to grab the `loggy` file to try to timestamp when each user was
active.

## Importing history

With these in hand, we can now actually import the `svn` history[^2] via:

```{code-block} bash
mkdir eon_svn && cd eon_svn
git svn init https://theory.cm.utexas.edu/svn/eon --no-metadata --prefix ""
git config svn.authorsfile ../authors.txt
git svn fetch
```

This will take much longer, since `fetch`, unlike `checkout` downloads all
changes. Cross-reference with the "fuller" variant from here:

```{code-block} bash
svn log | awk -F '|' '/^r/ {sub("^ ", "", $2); 
sub(" $", "", $2); print $2" = "$2" <"$2">"}' | sort -u > authors-transform.txt
# Includes things like
(no author) = (no author) <na@NODATA>
```

By project convention, if the command fails with an author not defined use
`username = username <username@NODATA>`, then retry the `git svn fetch`.

Also note that it is by far preferably to use Github handles for email
addresses, `<haozeke@users.noreply.github.com>` to reduce leaking personally
identifiable information.

## Merging upstream

Since the SVN source has no tags or branches, there are no cleanup tasks
required.

```{note}
`eON` used to use CVS for managing history, which was not imported when the switch to SVN took place, so rememeber to rebase a dummy commit co-authored by contributors before 2010!
```

Since we already have `svn` as a branch and even tags, it is a bit of a pain to
take updates. Essentially we will have to rebase onto the "new" commits.

```{code-block} bash
git branch -d main # extraneous
git checkout -b svn_new
git remote add origin git@github.com:TheochemUI/eongit
git pull
# Add in the commits from the older svn branch
# Just the gitignore bits
git cherry-pick b8794c316882b09e6421a3ba9f707a893accb636
git cherry-pick 35ea1448ee8707e06d4b314fbdfb0fabe859a206
```

Now we are ready to cut and push a new tag.

```{code-block} bash
git checkout svn_new
git checkout -b main_tk3
git reset --hard $LAST_TAG_COMMIT # July 1 this is 0e263123
# from to
git rebase --ignore-whitespace --rebase-merges --onto main_tk3 b8794c3^ 5505ae2
# After fixing many changes 594826a
git branch main_fin
git push -u origin main_fin
# Adding back the newer commits after the reset
git rebase --ignore-whitespace --rebase-merges --onto main_fin 0e26312^ c267e1c
git push origin HEAD:main_fin
```

Note that the `ignore-whitespace` option might be problematic for Python
changes, but for C++ only changesets it should be fine.

Documentation only commits are ported separately, and skipped during the rebase, since these will be in files "deleted by us".

```{code-block} bash
# Open next both merged
emacsclient -n $(git diff --name-only --diff-filter=U | head -n 1)
```


```{note}
Remember to diff the folders (e.g. via `meld`) after this operation, and confirm only expected changes are present.
```

[^1]: Atlassian has [pretty decent documentation and helpers](https://www.atlassian.com/git/tutorials/svn-to-git-prepping-your-team-migration), but, Java is a pain.
[^2]: From this [SO answer](https://stackoverflow.com/a/79188/1895378)
