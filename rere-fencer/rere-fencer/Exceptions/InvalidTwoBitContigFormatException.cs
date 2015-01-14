using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace rere_fencer.Exceptions
{
    public class InvalidTwoBitContigFormatException : Exception
    {
        public InvalidTwoBitContigFormatException(string contig, string message) 
            : base(string.Format("Contig {0} is not a valid TwoBit contig.\n{1}", contig, message)) { }
    }
}
