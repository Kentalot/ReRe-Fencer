using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.DataStructures.Contigs;
using BioinformaticUtils.DataStructures.Nucleotides;

namespace rere_fencer.Processors
{
    public class ParallelFastaWriter : IGenomeWriter
    {
        private readonly IEnumerable<IGenomeContig> _genomeContigs;
        private readonly FileInfo _outputFile;
        private const int MaxNucsPerLine = 50;

        public ParallelFastaWriter(IEnumerable<IGenomeContig> genomeContigs, FileInfo outputFilesLookup)
        {
            _genomeContigs = genomeContigs;
            _outputFile = outputFilesLookup;
        }

        public void WriteAllContigs()
        {
            var outputfile = _outputFile.FullName;
            if (outputfile.EndsWith(".fa"))
                outputfile = outputfile.Substring(0, outputfile.Length - 3);
            Parallel.ForEach(_genomeContigs, contig =>
            {
                var charBuffer = new char[MaxNucsPerLine];
                using (var writer = new StreamWriter(new FileStream(outputfile + "." + contig.Name + ".fa", FileMode.Create, FileAccess.Write)))
                {
                    writer.WriteLine(">" + contig.Name);
                    var i = 0;
                    foreach (var nuc in contig.GetWholeSequence())
                    {
                        if (i == MaxNucsPerLine)
                        {
                            writer.WriteLine(new string(charBuffer));
                            charBuffer = new char[MaxNucsPerLine];
                            i = 0;
                        }
                        charBuffer[i++] = nuc;
                    }
                    writer.WriteLine(new string(charBuffer.Take(i).ToArray()));
                }
            });
        }
    }
}
