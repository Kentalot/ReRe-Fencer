using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BioinformaticUtils.DataStructures.Contigs;

namespace rere_fencer.Processors
{
    public interface IRRFResolver
    {
        void Resolve(IEnumerable<IGenomeContig> genomeContigs, IGenomeWriter genomeWriter);
    }
}
