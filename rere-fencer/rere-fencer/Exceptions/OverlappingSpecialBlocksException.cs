using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace rere_fencer.Exceptions
{
    public class OverlappingSpecialBlocksException : Exception
    {
        public OverlappingSpecialBlocksException(string contig, string message) 
            : base(string.Format("Contig {0} has overlapping blocks.\n{1}", contig, message)) { }
    }
}
