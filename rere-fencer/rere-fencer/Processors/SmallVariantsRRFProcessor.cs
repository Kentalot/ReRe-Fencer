using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.DataStructures.Contigs;
using BioinformaticUtils.DataStructures.Nucleotides;
using BioinformaticUtils.Vcf;
using BioinformaticUtils.Vcf.VariantInfoClasses;

namespace rere_fencer.Processors
{
    public class SmallVariantsRRFProcessor : IRRFProcessor
    {
        private static readonly char missingChar = VcfVariant.MissingValueString.First();
        public IEnumerable<DnaNucleotide> Process(IGenomeContig genomeContig, IEnumerable<IVcfVariant> vcfVariants)
        {
            var posOffset = genomeContig.IsZeroBasedCoordinates ? 0U : 1U;
            uint contigPosition = 0U;

            foreach (var variantInfo in vcfVariants.Select(v => v.CreateVariantInfo()).Where(v => v.IsPassFilter && !v.IsRefCall))
            {
                var variant = variantInfo.Variant;
                if (contigPosition + posOffset > variant.Position) continue; // overlapping variants are ignored, take the first one only.

                //check for any hemi or hom.
                var hemiOrHom = 0;
                for (; hemiOrHom < variant.Samples.Count; hemiOrHom++)
                    if (variantInfo.SampleInfo.IsHemi(hemiOrHom) || variantInfo.SampleInfo.IsHom(hemiOrHom))
                        break;
                if (hemiOrHom == variant.Samples.Count) continue; // if none, then go to next variant.
                var altIndex = int.Parse(variant.Samples[hemiOrHom][VcfVariant.GenotypeKey][0].ToString(CultureInfo.CurrentCulture)) - 1;
                if (variant.Alts[altIndex].Any(c => !VcfVariant.ValidAltNucleotides.Contains<DnaNucleotide>(c)))
                {
                    Console.WriteLine("Skipping variant with unrecognized Alt: " + variant);
                    continue;
                }

                uint end = 0;
                if (variant.Position != 1) // edge case if it's 1, then just take the alt alleles.
                {
                    end = variant.Position + posOffset - 2; // take sequence up to base before first ref base
                    if (end >= contigPosition)
                        // could happen where there's a snp and then a snp that ends up making end < contigPosition
                        foreach (var nuc in genomeContig.GetSequence(contigPosition, end))
                            yield return nuc;
                    end++; // add 1 to end so that when contigPosition is calculated, it ends up 1 more than the length of ref.
                }

                contigPosition = end + (uint)variant.Ref.Length; // no need for posOffset since that was accounted for in end.

                if (variant.Alts[altIndex].Equals(VcfVariant.MissingValueString)) continue;
                foreach (var nuc in variant.Alts[altIndex])
                    yield return nuc;
            }

            foreach (var nuc in genomeContig.GetSequence(contigPosition, (uint)genomeContig.Length + posOffset - 1))
                yield return nuc;
        }
    }
}
