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
            uint contigPosition = 0U;

            foreach (var variantInfo in vcfVariants.Select(v => v.CreateVariantInfo()))
            {
                var variant = variantInfo.Variant;
                if (contigPosition + posOffset > variant.Position) continue; // overlapping variants are ignored, take the first one only.

                //check for any hemi or hom.
                var hemiOrHom = 0;
                for (; hemiOrHom < variant.Samples.Count; hemiOrHom++)
                    if (variantInfo.SampleInfo.IsHemi(hemiOrHom) || variantInfo.SampleInfo.IsHom(hemiOrHom))
                        break;
                if (hemiOrHom == variant.Samples.Count) continue; // if none, then go to next variant.

                var end = variant.Position + posOffset - 2;// take sequence up to base before first ref base
                foreach (var nuc in genomeContig.GetSequence(contigPosition, end))
                {
                    yield return nuc;
                    //contigPosition[0]++;
                }
                contigPosition = end + 1 + (uint)variant.Ref.Length; // no need for posOffset since that was accounted for in end.

                var altIndex = int.Parse(variant.Samples[hemiOrHom][VcfVariant.GenotypeKey][0].ToString(CultureInfo.CurrentCulture)) - 1;
                foreach (var nuc in variant.Alt[altIndex])
                    yield return nuc;  
            }

            foreach (var nuc in genomeContig.GetSequence(contigPosition, (uint)genomeContig.Length + posOffset - 1))
                yield return nuc;
        }
    }
}
