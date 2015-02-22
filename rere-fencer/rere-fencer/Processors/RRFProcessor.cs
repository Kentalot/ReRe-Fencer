using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.DataStructures.Nucleotides;
using BioinformaticUtils.GenomeTools;
using BioinformaticUtils.Vcf;
using BioinformaticUtils.Vcf.VariantInfoClasses;

namespace rere_fencer.Processors
{
    public class RRFProcessor : IRRFProcessor
    {

        public IEnumerable<DnaNucleotide> Process(IGenomeReader genomeReader, IVcfReader vcfReader)
        {
            var posOffset = genomeReader.IsZeroBasedCoordinates ? 0U : 1U;
            uint[] position = {0U};
            string currentcontig = null;
            foreach (var variantInfo in vcfReader.ReadVcfVariants().Select(v => v.CreateVariantInfo()))
            {
                var hemiOrHom = -1;
                if (currentcontig != variantInfo.Variant.Chromosome.Name)
                    currentcontig = variantInfo.Variant.Chromosome.Name;
                for (var i = 0; i < variantInfo.Variant.Samples.Count; i++)
                    if (variantInfo.SampleInfo.IsHemi(i) || variantInfo.SampleInfo.IsHom(i))
                    {
                        hemiOrHom = i;
                        break;
                    }
                if (hemiOrHom < 0) continue;

                var end = variantInfo.Variant.Position + (uint) variantInfo.Variant.Ref.Length - 1;
                var seq = genomeReader.Contigs.GetSequence(currentcontig, position[0], end);
                var info = variantInfo;
                foreach (var nuc in seq.TakeWhile(nuc => position[0] + posOffset != info.Variant.Position))
                {
                    yield return nuc;
                    position[0]++;
                }
                position[0] = end + 1;

                var altIndex = int.Parse(variantInfo.Variant.Samples[hemiOrHom][VcfVariant.GenotypeKey][0].ToString(CultureInfo.CurrentCulture));
                foreach (var nuc in info.Variant.Alt[altIndex])
                    yield return nuc;
            }
        }
    }
}
