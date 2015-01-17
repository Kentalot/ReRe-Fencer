using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace rere_fencer.Exceptions
{
    public class OutOfContigRangeException : Exception
    {
        public OutOfContigRangeException(string contig, uint start, uint end, string message) 
            : base(string.Format("Tried to access a sequence ({2}, {3}) completely outside of the Contig: {0} .\n{1}", 
                                 contig, message, start, end)) { }
    }

    public class OutOfSubSequenceRangeException : ArgumentOutOfRangeException
    {
        public OutOfSubSequenceRangeException(uint start, uint end, string type) 
            : base(string.Format("Tried to access a block of subsequence great than 2,147,483,644 in length!\n" +
                "SubSequence Start: {0}; End: {1}; Type: {2}", start, end, type)) { }
    }
}
