using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.DataStructures.Contigs;
using BioinformaticUtils.GenomeTools;

namespace rere_fencer.Processors
{
    public interface IGenomeWriter
    {
        void TryWriteContigToStream(IGenomeContig genomeContig);
    }
}
