using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.DataStructures.Contigs;

namespace rere_fencer.Processors
{
    public class GenomeWriter : TextWriter, IGenomeWriter
    {
        public TextWriter WriteContig(IGenomeContig genomeContig)
        {
            throw new NotImplementedException();
        }

        public override Encoding Encoding
        {
            get { throw new NotImplementedException(); }
        }
    }
}
