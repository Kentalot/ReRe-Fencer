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
        IReadOnlyDictionary<IContig, IGenomeContig> Contigs { get; }
    }

    public interface IContig
    {
        string Name { get; }
        uint? Length { get; }
    }

    public sealed class Contig : IContig
    {
        public string Name { get; private set; }
        public uint? Length { get; private set; }

        public Contig(string name, uint? length)
        {
            Name = name;
            Length = length;
        }

        public Contig(string name) : this(name, null) { }

        public override bool Equals(Object obj)
        {
            var other = obj as Contig;
            return other != null && Equals(other);
        }

        public bool Equals(Contig other)
        {
            return string.Equals(Name, other.Name) && Length == other.Length;
        }

        public override int GetHashCode()
        {
            return (91) + Name.GetHashCode();
        }

        public static bool operator ==(Contig a, Contig b)
        {
            if (ReferenceEquals(a, b))
            {
                return true;
            }

            // If one is null, but not both, return false.
            if (((object)a == null) || ((object)b == null))
            {
                return false;
            }

            return a.Equals(b);
        }

        public static bool operator !=(Contig a, Contig b)
        {
            return !(a == b);
        }

        public override String ToString()
        {
            return Name;
        }
    }

    public interface IGenomeContig : IContig
    {
        bool ContainsNs { get; }
        bool ContainsMaskedSequences { get; }
        IEnumerable<char> GetSequence(uint start, uint end, ReadMode readMode = ReadMode.Normal);
        char GetNucleotideAt(uint position);
        IEnumerable<IContigRegion> MaskedRegions { get; }
        IEnumerable<IContigRegion> NRegions { get; }
        IEnumerable<IContigRegion> NormalRegions { get; }
        //uint? FirstPositionOfSequence(string sequence, uint offset = 0, bool strict = true);
        //uint? LastPositionOfSequence(string sequence, uint offset = 0, bool strict = true);
    }

    public enum ReadMode
    {
        Normal = 0, IgnoreMasks, SkipNs, IgnoreMasksSkipNs, SkipMasks, SkipMasksSkipNs
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
