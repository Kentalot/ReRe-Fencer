using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.GenomeTools;
using BioinformaticUtils.Vcf;

namespace rere_fencer.Processors
{
    public class RRFProcessor : IRRFProcessor
    {

        public void Process(IGenomeReader genomeReader, IVcfReader vcfReader)
        {
            var position = genomeReader.IsZeroBasedCoordinates ? 0U : 1U;
            foreach (var variantInfo in vcfReader.ReadVcfVariants().Select(v => v.CreateVariantInfo())
                .Where(v => v.SampleInfo.IsHom() || v.SampleInfo.IsHemi()))
            {
                var end = variantInfo.Variant.Position;
            }
        }
    }
}
