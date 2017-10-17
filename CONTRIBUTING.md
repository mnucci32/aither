# How to Contribute
Thank you for contributing to aither! Check out the aither [blog](http://aithercfd.com) for more information. Aither
uses the git flow branching method. The **master** branch is a permanent branch containing stable code with releases
coming off of this branch. The default branch is the **develop** branch. This is a permanent branch containing the most 
up-to-date code and is therefore not guaranteed to be working. Bug fixes and feature additions branch off of **develop** and 
are merged back in when complete. Aither uses continuous integration services [TravisCI](https://travis-ci.org/mnucci32/aither) 
and [Appveyor](https://ci.appveyor.com/project/mnucci32/aither/branch/develop) to test builds on Linux, macOS, and 
Windows. The code is built on each platform and a suite of tests are run to ensure that the residuals have not 
changed due to a commit. There are a few ways to get involved with this project including reporting or fixing a bug, 
requesting or adding a feature, and adding documentation.

### Reporting a Bug
Bugs should be reported via [Github Issues](https://github.com/mnucci32/aither/issues).

### Fixing a Bug
A list of known bugs can be found on the [Issues](https://github.com/mnucci32/aither/issues) page. A fix can be submitted 
by creating a new branch off of the **develop** branch with a descriptive name starting with **hotfix_**. When the bug has
been fixed and and all continuous integration tests are passing, a pull request can be submitted to merge the fix back into
**develop**.

### Requesting a Feature
New features should be requested via [Github Issues](https://github.com/mnucci32/aither/issues). Not all requests will be
accepted.

### Adding a Feature
A list of features to be added can be found on the [Issues](https://github.com/mnucci32/aither/issues) page. If you want
to add a feature not listed, first submit an issue requesting that it be added. New features should be implemented by 
creating a branch off of **develop** with a descriptive name starting with **feature_**. When the feature has been fixed 
and all continuous integreation tests are passing, a pull request can be submitted to merge the feature back into **develop**.
If it is a fairly significant addition, an additional test case should be added as well.

### Adding Documentation
Documentation can be added by editing the [Wiki](https://github.com/mnucci32/aither/wiki) page.

# Coding Style
* Aither conforms to the [Google Style Guide](https://google.github.io/styleguide/cppguide.html)
* Variable names should use camelCase starting with a lower case letter
* Function names should use CamelCase starting with an upper case letter
* Class member variables should be appended with a _
* Modern C++ features should be used when possible
