---
myst:
  html_meta:
    "description": "Information on obtaining older versions of EON from the Subversion (SVN) repository."
    "keywords": "EON SVN, Subversion, legacy code, version control"
---

# Subversion

```{deprecated} 2.0
These instructions are for obtaining older versions of `eOn`.
```

A read-only copy of the code can be anonymously checked out of the subversion
repository with the following command:

```{code-block} bash
svn co http://theory.cm.utexas.edu/svn/eon
```

If you have developer access to the UT Austin hosted SVN instance, you can
checkout the code using the command:

```{code-block} bash
svn co svn+ssh://username@theory.cm.utexas.edu/svn/eon
```

This will fetch a copy of the latest code to a local directory named eon.

For more information on subversion read its [online documentation](http://svnbook.red-bean.com/en/1.5/index.html).

The `eOn` code can also be download directly in `tar` format [from here](http://theory.cm.utexas.edu/code/eon.tgz).
