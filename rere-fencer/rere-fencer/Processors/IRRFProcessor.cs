using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.GenomeTools;
using BioinformaticUtils.Vcf;

namespace rere_fencer.Processors
{
    public interface IRRFProcessor
    {
        void Process(IGenomeReader genomeReader, IVcfReader vcfReader);
    }
}
