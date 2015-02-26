using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.DataStructures.Contigs;
using BioinformaticUtils.DataStructures.Nucleotides;
using BioinformaticUtils.GenomeTools;
using BioinformaticUtils.Vcf;
using BioinformaticUtils.Vcf.VariantInfoClasses;

namespace rere_fencer.Processors
{
    public class SmallVariantsRRFProcessor : IRRFProcessor
    {
        public IEnumerable<DnaNucleotide> Process(IGenomeContig genomeContig, IEnumerable<IVcfVariant> vcfVariants)
        {
            var posOffset = genomeContig.IsZeroBasedCoordinates ? 0U : 1U;
            uint[] contigPosition = {0U};

            foreach (var variantInfo in vcfVariants.Select(v => v.CreateVariantInfo()))
            {
                if (contigPosition[0] + posOffset > variantInfo.Variant.Position) continue; // overlapping variants are ignored, take the first one only.

                //check for any hemi or hom.
                var hemiOrHom = 0;
                for (; hemiOrHom < variantInfo.Variant.Samples.Count; hemiOrHom++)
                    if (variantInfo.SampleInfo.IsHemi(hemiOrHom) || variantInfo.SampleInfo.IsHom(hemiOrHom))
                        break;
                if (hemiOrHom == variantInfo.Variant.Samples.Count) continue; // if none, then go to next variant.

                var end = variantInfo.Variant.Position + (uint)variantInfo.Variant.Ref.Length - 1U;
                var seq = genomeContig.GetSequence(contigPosition[0], end);
                var info = variantInfo;
                foreach (var nuc in seq) // assume that the vcf ref sequence matches with the genome ref sequence.
                {
                    yield return nuc;
                    contigPosition[0]++;
                }
                contigPosition[0] = end + 1;

                var altIndex = int.Parse(variantInfo.Variant.Samples[hemiOrHom][VcfVariant.GenotypeKey][0].ToString(CultureInfo.CurrentCulture));
                foreach (var nuc in info.Variant.Alt[altIndex])
                    yield return nuc;  
            }

            foreach (var nuc in genomeContig.GetSequence(contigPosition[0], (uint)genomeContig.Length + posOffset - 1))
                yield return nuc;
        }
    }
}
