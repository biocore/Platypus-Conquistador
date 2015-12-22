# Platypus-Conquistador changelog

Platypus Conquistador 0.9.0-dev (changes since Platypus Conquistador 0.9.0 go here)
-----------------------------------------------------------------------------------

- Changed parse.parse_m9 so it doesn't load the full file into memory
- Memory requirements was improved by moving the labels and counts of hits to the databases
to simple files with all the labels of the hits. Basically, now those files have all the
labels of the hit sequences. For a summary, you could simply do: `cat your_file | sort | uniq -c`.

Version 0.9.0 (2015-04-26)
--------------------------

Initial beta release.
