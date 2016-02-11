# Prerequisites

This analysis depends on the following software.

- [vmatch][vmatch]
- [GenHub][genhub], which in turn depends on:
    - [AEGeAn][agn]
    - [GenomeTools][gt]

The installation instructions for vmatch, AEGeAn, and GenomeTools should be sufficient.
However, the default GenHub installation with pip will not work currently.
The recipes for TAIR6, OGS1.0, and OGS3.2 are only available from the `robust` branch of the GitHub repository.
Install GenHub as follows.

```bash
git clone https://github.com/standage/genhub.git
cd genhub
git checkout robust
python setup.py install
```

GenHub's [installation guide][genhub-install] *does* include some helpful troubleshooting tips if you run into issues.


[vmatch]: http://www.vmatch.de/
[genhub]: https://github.com/standage/genhub
[agn]: https://brendelgroup.github.io/AEGeAn
[gt]: https://genometools.org
[genhub-install]: https://github.com/standage/genhub/blob/master/docs/INSTALL.md
