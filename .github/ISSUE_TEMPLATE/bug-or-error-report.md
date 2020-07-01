---
name: Bug or Error report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug or error is.

**Error message**
Paste the complete log message, include the main task log and failed subtask log.
The main task log is usually located  in your working directory and is named `pidXXX.log.info` and the main task log will tell you the failed subtask log in the last few lines, such as:
<pre>
[ERROR] 2020-07-01 11:06:57,184 cns_align failed: please check the following logs:
[ERROR] 2020-07-01 11:06:57,185 <ins>~/NextDenovo/test_data/01_rundir/02.cns_align/02.cns_align.sh.work/cns_align0/nextDenovo.sh.e</ins>
</pre>

**Operating system**
Which operating system and version are you using?
You can use the command `lsb_release -a` to get it.

**GCC**
What version of GCC are you using?
You can use the command `gcc -v` to get it.

**Python**
What version of Python are you using?
You can use the command `python --version` to get it.

**NextDenovo**
What version of NextDenovo are you using?
You can use the command `nextDenovo -v` to get it.

**To Reproduce** (Optional)
Steps to reproduce the behavior. Providing a minimal test dataset on which we can reproduce the behavior will generally lead to quicker turnaround time!

**Additional context** (Optional)
Add any other context about the problem here.
