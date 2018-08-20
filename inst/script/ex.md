This file describes the 'ex' directory

'ex' is a toy example for the vignette and users to play around with.
The samples PCR1, PCR2, and PCR3 are generated using synthetic datasets.

The dataset used in the above examples
was obtained from a combination of synthetic sample datasets generated
using [MiXCR](https://github.com/milaboratory/mixcr)'s program 
[here](http://files.milaboratory.com/mixcr/paper/mixcr-test-1.2-SNAPSHOT.jar). 
Firstly, three distinct samples were generated, each simulated with the following parameters
in `MiXCR`:

```{r mixr-test-params, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'}
tabl <- "
| Parameter     | sample 1      | sample 2      | sample 3      |
|---------------|---------------|---------------|---------------|
| reads         |  10000        |  10000        | 10000         |
| clones        |  5000         |  5000         | 2000          |
| seed          |  4228         |  2428         | 2842          |
| conf          | MiSeq-300-300 | MiSeq-300-300 | MiSeq-300-300 |
| loci          |  IGH          | IGH           | IGH           |
| species       | hsa           | hsa           | hsa           |
"
cat(tabl)
```

Following that, an arbitrary number of sequences were randomly
drawn from each of the 3 samples and randomly amplified.
This process was repeated 3 times, resulting in a final repertoire
of 3 samples, named PCR1, PCR2, and PCR3.

Finally, these 3 samples were analyzed by sabseqPy. The command used to
analyze these samples are as follows:

```bash
abseq -y params.yml
```

where the contents of `params.yml` is:

```yaml
# params.yml
defaults:
    bitscore: 300
    sstart: 1-3
    alignlen: 250
    outdir: ex
    task: all
    threads: 1
---
file1: PCR1.fasta
name: PCR1
---
file1: PCR2.fasta
name: PCR2
---
file1: PCR3.fasta
name: PCR3
```

abseqPy's analysis output on these 3 samples are contained
within the dataset described in this vignette.

