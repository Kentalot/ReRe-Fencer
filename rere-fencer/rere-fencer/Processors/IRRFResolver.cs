using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace rere_fencer.Processors
{
    public interface IRRFResolver
    {
        void Resolve(IRRFProcessor rrfProcessor, IGenomeWriter genomeWriter);
    }
}
