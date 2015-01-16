using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace rere_fencer.Input
{
    public interface IGenomeReader
    {
        bool IsZeroBasedCoordinates { get; }
        bool SupportsMasking { get; }
        bool SupportsNs { get; }
        bool SupportsIupacAmbiguityCodes { get; }
        SortedList<string, IGenomeContig> Contigs { get; }
    }

    public interface IGenomeContig
    {
        string Name { get; }
        uint Length { get; }
        bool ContainsNs { get; }
        bool ContainsMaskedSequences { get; }
        string GetSequenceAt(uint start, uint end, bool ignoreNs = false);
        char GetNucleotideAt(uint position);
        uint? FirstPositionOfSequence(string sequence, uint offset = 0, bool strict = true);
        uint? LastPositionOfSequence(string sequence, uint offset = 0, bool strict = true);
    }
}
