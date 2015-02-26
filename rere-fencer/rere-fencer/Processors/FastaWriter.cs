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
    public class FastaWriter : IGenomeWriter
    {
        private readonly IEnumerable<IGenomeContig> _genomeContigs;
        private readonly FileInfo _outputFile;
        private const int MaxNucsPerLine = 50;

        public FastaWriter(IEnumerable<IGenomeContig> genomeContigs, FileInfo outputFilesLookup)
        {
            _genomeContigs = genomeContigs;
            _outputFile = outputFilesLookup;
        }

        public void WriteAllContigs()
        {
            var charBuffer = new char[MaxNucsPerLine];
            using (var writer = new StreamWriter(new FileStream(_outputFile.FullName, FileMode.Create, FileAccess.Write)))
            {
                foreach (var contig in _genomeContigs)
                {
                    writer.WriteLine("> " + contig.Name);
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
            }
        }
    }
}
