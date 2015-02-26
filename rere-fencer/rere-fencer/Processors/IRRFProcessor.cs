using System;
using System.Collections.Generic;
using BioinformaticUtils.DataStructures.Contigs;
using BioinformaticUtils.DataStructures.Nucleotides;
using BioinformaticUtils.GenomeTools;
using BioinformaticUtils.Vcf;

namespace rere_fencer.Processors
{
    public interface IRRFProcessor
    {
        IEnumerable<DnaNucleotide> Process(IGenomeContig genomeContig, IEnumerable<IVcfVariant> vcfVariants);
    }
}
