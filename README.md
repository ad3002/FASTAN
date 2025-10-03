# FasTAN: A Fast Tandem Repeat Finder  
<font size ="4">**_Author:  Gene Myers_**<br>
**_First:   Sept 30, 2025_**<br>

- [FasTAN](#FasTAN) 

## Overview

**FasTAN** blah blah blah

<a name="FasTAN"></a>

## FasTAN Reference

```
FasTAN <source:path>[<precursor] <target:path>[.1aln]
                  
    <precursor> = .1gdb | <fa_extn> | <1_extn>
    
    <fa_extn> = (.fa|.fna|.fasta)[.gz]
    <1_extn>  = any valid 1-code sequence file type
```

FasTAN takes a FASTA or 1GDB as input and outputs the off-diagonal alignment delimiting
tandemly repeats and their period as a .1aln file.  A prototype, doesn't even make
cleanly.