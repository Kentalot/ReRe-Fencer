using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace rere_fencer.Exceptions
{
    public class InvalidTwoBitFileException : Exception
    {
        public InvalidTwoBitFileException(string file, string message) 
            : base(file + " file is not a valid TwoBit file.\n" + message) { }
    }
}
