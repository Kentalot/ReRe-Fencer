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
        IReadOnlyList<IGenomeContig> Contigs { get; }
    }

    public interface IGenomeContig
    {
        string Name { get; }
        uint Length { get; }
        bool ContainsNs { get; }
        bool ContainsMaskedSequences { get; }
        IEnumerable<char> GetSequence(uint start, uint end, bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false);
        char GetNucleotideAt(uint position);
        IEnumerable<IContigRegion> MaskedRegions { get; }
        IEnumerable<IContigRegion> NRegions { get; }
        IEnumerable<IContigRegion> NormalRegions { get; }
        //uint? FirstPositionOfSequence(string sequence, uint offset = 0, bool strict = true);
        //uint? LastPositionOfSequence(string sequence, uint offset = 0, bool strict = true);
    }

    public interface IContigRegion
    {
        uint Start { get; }
        uint End { get; }
    }

    internal class ContigRegion : IContigRegion
    {
        private readonly Tuple<uint, uint> _regionInfo;
        public uint Start { get { return _regionInfo.Item1; } }
        public uint End { get { return _regionInfo.Item2; } }
        public ContigRegion(uint start, uint end) { _regionInfo = Tuple.Create(start, end); }
    }
}
