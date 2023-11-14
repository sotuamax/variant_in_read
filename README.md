# variant_in_read
check variants in alignment reads compared to reference fasta 

### input files:
Required: 
1. reference fasta file (REF)
2. BAM alignment file (BAM)
3. output name (OUT) <br>
Optional: <br>
BED file <br>

To enable multiple cores processing, use:
```
mpiexec -n 8 fidelity_check.py REF BAM -o OUT 
```
