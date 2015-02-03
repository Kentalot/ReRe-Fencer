using System;
using System.CodeDom;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Globalization;
using System.IO;
using System.IO.MemoryMappedFiles;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Schema;
using rere_fencer.Exceptions;

namespace rere_fencer.Input
{

    public class TwoBitGenomeReader : IGenomeReader, IDisposable
    {
        private class TwoBitGenomeContig : IGenomeContig, IDisposable
        {
            public string Name { get; private set; }
            public uint Length { get; private set; }
            public uint Reserved { get; private set; }
            public bool ContainsNs { get { return NRegions.Any(); } }
            public bool ContainsMaskedSequences { get { return MaskedRegions.Any(); } }
            public IEnumerable<IContigRegion> NRegions 
                { get { return _nBlockStarts.Zip(_nBlockEnds, (x, y) => new ContigRegion(x, y)); } }
            public IEnumerable<IContigRegion> MaskedRegions
            { get { return _maskBlockStarts.Zip(_maskBlockEnds, (x, y) => new ContigRegion(x, y)); } }

            private IList<IContigRegion> _normalRegions; 
            public IEnumerable<IContigRegion> NormalRegions
            { 
                get { return GetNormalRegions(); } 
            }

            //private readonly List<ITwoBitGenomeSubcontig> _subcontigs;
            private readonly MemoryMappedViewAccessor _contigAccessor;
            private readonly uint[] _nBlockStarts;
            private readonly uint[] _nBlockEnds;
            private readonly uint[] _maskBlockStarts;
            private readonly uint[] _maskBlockEnds;
            
            public TwoBitGenomeContig(string name, uint offset, uint size)//uint length, List<ITwoBitGenomeSubcontig> subcontigs)
            {
                Name = name;
                var bytePosition = 4U;
                using (var tempViewer = _genomeFile.CreateViewAccessor(offset, size, MemoryMappedFileAccess.Read))
                {
                    Length = ReadUint(tempViewer); // read first 4 bytes.
                    var numBlocks = ReadUint(tempViewer, bytePosition); // read 4 bytes in.
                    bytePosition += 4;
                    _nBlockStarts = new uint[numBlocks];
                    _nBlockEnds = new uint[numBlocks];
                    bytePosition = GetBlockStartsAndEnds(tempViewer, bytePosition, _nBlockStarts, _nBlockEnds);
                    numBlocks = ReadUint(tempViewer, bytePosition);
                    bytePosition += 4;
                    _maskBlockStarts = new uint[numBlocks];
                    _maskBlockEnds = new uint[numBlocks];
                    bytePosition = GetBlockStartsAndEnds(tempViewer, bytePosition, _maskBlockStarts, _maskBlockEnds);
                    Reserved = ReadUint(tempViewer, bytePosition);
                    bytePosition += 4;
                }
                _contigAccessor = _genomeFile.CreateViewAccessor(bytePosition + offset, size == 0 ? 0 : size - bytePosition,
                    MemoryMappedFileAccess.Read);
                //Length = length;
                //_subcontigs = subcontigs;
            }

            public IEnumerable<char> GetSequence(uint start, uint end, ReadMode readMode = ReadMode.Normal)
            {
                if (end >= Length)
                    throw new OutOfContigRangeException(Name, start, end,
                        string.Format("Inside {0}'s GetSequence method.", GetType().Name));
                if (end < start) throw new ArgumentException("Sequence end was less than start, what's up with that?");
                switch (readMode)
                {
                    case ReadMode.SkipMasks:
                        return GetSequenceWithOptions(start, end, skipMasks: true);
                    case ReadMode.IgnoreMasks:
                        return GetSequenceWithOptions(start, end, true);
                    case ReadMode.SkipNs:
                        return GetSequenceWithOptions(start, end, skipNs: true);
                    case ReadMode.SkipMasksSkipNs:
                        return GetSequenceWithOptions(start, end, false, true, true);
                    case ReadMode.IgnoreMasksSkipNs:
                        return GetSequenceWithOptions(start, end, true, true);
                    default:
                        return GetSequenceWithOptions(start, end);
                }
            }

            private IEnumerable<char> GetSequenceWithOptions(uint start, uint end, bool ignoreMasks = false, bool skipNs = false, bool skipMasks = false)
            {
                var length = end - start + 1;
                uint charPosition = 0, bytePosition = start/NucsPerByte;
                var leftStartIndex = start % NucsPerByte;
                var numBytes = (int) Math.Ceiling((length + leftStartIndex)/(decimal)NucsPerByte);
                var overlappingNRange = GetOverlappingBlocksRange(start, end, true);
                var overlappingMaskRange = ignoreMasks ? null : GetOverlappingBlocksRange(start, end, false);
                bool noMoreNs = false, noMoreMasks = false;
                int nIndex, maskIndex;
                if (overlappingNRange == null)
                {
                    noMoreNs = true;
                    nIndex = -1;
                }
                else nIndex = overlappingNRange.Item1;
                if (overlappingMaskRange == null)
                {
                    noMoreMasks = true;
                    maskIndex = -1;
                }
                else maskIndex = overlappingMaskRange.Item1;
                var position = start;
                var inNBlock = IsOverlapping(position, _nBlockStarts, _nBlockEnds, nIndex);
                var inMaskBlock = IsOverlapping(position, _maskBlockStarts, _maskBlockEnds, maskIndex);
                for (var j = 0; j < numBytes; j++)
                {
                    var tetNuc = ByteToTetraNucleotide(_contigAccessor.ReadByte(bytePosition));
                    for (var i = j == 0 ? (int) leftStartIndex : 0; i < NucsPerByte; i++)
                    {
                        if (!noMoreNs && inNBlock)
                        {
                            if (!skipNs) yield return 'N';
                            if (++charPosition == length) yield break;
                            position++;
                            inNBlock = IsOverlapping(position, _nBlockStarts, _nBlockEnds, nIndex);
                            if (inNBlock) continue;
                            if (++nIndex == overlappingNRange.Item2 + 1) noMoreNs = true; // already at the end of the nBlocks
                            else inNBlock = IsOverlapping(position, _nBlockStarts, _nBlockEnds, nIndex);
                        }
                        else if (!noMoreMasks && inMaskBlock)
                        {
                            if (!skipMasks) yield return char.ToLower(tetNuc[i], CultureInfo.CurrentCulture);
                            if (++charPosition == length) yield break;
                            position++;
                            inMaskBlock = IsOverlapping(position, _maskBlockStarts, _maskBlockEnds, maskIndex);
                            if (inMaskBlock) continue;
                            if (++maskIndex == overlappingMaskRange.Item2 + 1) noMoreMasks = true; // already at the end of maskBlocks
                            else inMaskBlock = IsOverlapping(position, _maskBlockStarts, _maskBlockEnds, maskIndex);
                        }
                        else
                        {
                            yield return tetNuc[i];
                            if (++charPosition == length) yield break;
                            position++;
                            if (!noMoreNs) inNBlock = IsOverlapping(position, _nBlockStarts, _nBlockEnds, nIndex);
                            else if (!noMoreMasks) inMaskBlock = IsOverlapping(position, _maskBlockStarts, _maskBlockEnds, maskIndex); 
                        }
                    }
                    bytePosition++;
                }
            }

            private bool IsOverlapping(uint position, uint[] blockStarts, uint[] blockEnds, int index)
            {
                return index >= 0 && position >= blockStarts[index] && position <= blockEnds[index];
            }

            private Tuple<int, int> GetOverlappingBlocksRange(uint start, uint end, bool searchNs)
            {
                var starts = searchNs ? _nBlockStarts : _maskBlockStarts;
                var ends = searchNs ? _nBlockEnds : _maskBlockEnds;
                if (!starts.Any()) return null;
                var firstPossibleOverlap = BinarySearchForFirstOverlappingBlock(starts, ends, start, 0, starts.Length);
                if (firstPossibleOverlap < 0 || end < starts[firstPossibleOverlap]) // no overlaps possible.
                    return null;
                var lastPossibleOverlap = BinarySearchForLastOverlappingBlock(starts, ends, end, firstPossibleOverlap, starts.Length);
                if (lastPossibleOverlap < 0) return Tuple.Create(firstPossibleOverlap, firstPossibleOverlap);
                return Tuple.Create(firstPossibleOverlap, lastPossibleOverlap);
            }

            private int BinarySearchForLastOverlappingBlock(uint[] starts, uint[] ends, uint position, int begin, int end)
            {
                while (true)
                {
                    if (begin == end) return begin;
                    var middle = (begin + end)/2;
                    if (position < starts[middle]) // position is before the middle start
                    {
                        if (middle == 0) return -1;
                        end = middle - 1;
                        continue;
                    }
                    if (position > ends[middle])
                    {
                        if (middle == starts.Length - 1 // last one must be the last possible overlap.
                            || position < starts[middle + 1]) // this one must be the last possible overlap.
                            return middle;
                        begin = middle + 1;
                        continue;
                    }
                    return middle;
                    break;
                }
            }

            private int BinarySearchForFirstOverlappingBlock(uint[] starts, uint[] ends, uint position, int begin, int end)
            {
                while (true)
                {
                    var middle = (begin + end)/2;
                    if (position > ends[middle]) // position is after the middle end
                    {
                        if (middle + 1 == starts.Length) return -1; // cannot have any overlap.
                        begin = middle + 1;
                        continue;
                    }
                    if (position < starts[middle]) // position is before the middle start
                    {
                        if (middle == 0 // first one must be the first possible overlap.
                            || position > ends[middle - 1]) // this one must be the first possible overlap
                            return middle;
                        end = middle - 1;
                        continue;
                    }
                    return middle;
                    break;
                }
            }

            public char GetNucleotideAt(uint position)
            {
                return GetSequence(position, position).First();
            }

            private uint GetBlockStartsAndEnds(MemoryMappedViewAccessor tempViewer, uint start, uint[] blockStarts, uint[] blockEnds)
            {
                const uint uintByteCount = 4;
                var blockCount = blockStarts.Length;
                
                Parallel.For(0, blockCount, i =>
                {
                    blockStarts[i] = ReadUint(tempViewer, (uint) (i*uintByteCount + start));
                    blockEnds[i] = blockStarts[i] + ReadUint(tempViewer, (uint) (i*uintByteCount + blockCount*uintByteCount + start)) - 1;
                });

                return (uint) (start + blockCount*uintByteCount*2);
            }

            private IEnumerable<IContigRegion> GetNormalRegions()
            {
                if (_normalRegions == null)
                {
                    _normalRegions = new List<IContigRegion>();
                    uint nextBlock = 0, currentSequencePosition = 0;
                    for (int currentNIndex = 0, currentMaskIndex = 0;
                        currentNIndex < _nBlockStarts.Length || currentMaskIndex < _maskBlockStarts.Length; )
                    {
                        if (currentNIndex < _nBlockStarts.Length &&
                            currentSequencePosition == _nBlockStarts[currentNIndex])
                        {
                            nextBlock = _nBlockEnds[currentNIndex++];
                        }
                        else if (currentMaskIndex < _maskBlockStarts.Length &&
                                 currentSequencePosition == _maskBlockStarts[currentMaskIndex])
                        {
                            nextBlock = _maskBlockEnds[currentMaskIndex++];
                        }
                        else
                        {
                            var nBlockKey = currentNIndex == _nBlockStarts.Length
                                ? Length - 1
                                : _nBlockStarts[currentNIndex];
                            var maskBlockKey = currentMaskIndex == _maskBlockStarts.Length
                                ? Length - 1
                                : _maskBlockStarts[currentMaskIndex];
                            if (currentSequencePosition < nBlockKey && currentSequencePosition < maskBlockKey)
                            {
                                nextBlock = maskBlockKey < nBlockKey ? maskBlockKey : nBlockKey;
                                _normalRegions.Add(new ContigRegion(currentSequencePosition, nextBlock - 1));
                            }
                        }
                        currentSequencePosition = nextBlock;
                    }

                    var lastOne = _normalRegions.LastOrDefault();
                    if (lastOne == null || lastOne.End != Length - 1)
                        // means that there's one more normal block left to make
                        _normalRegions.Add(new ContigRegion(currentSequencePosition, Length - 1));
                }
                return _normalRegions;
            }
            
            public void Dispose()
            {
                _contigAccessor.Dispose();
            }
        }

        public bool IsZeroBasedCoordinates { get { return true; } }
        public bool SupportsMasking { get { return true; } }
        public bool SupportsNs { get { return true; } }
        public bool SupportsIupacAmbiguityCodes { get { return false; } }

        private readonly IGenomeContig[] _contigs;
        public IReadOnlyList<IGenomeContig> Contigs { get { return Array.AsReadOnly(_contigs); } }

        private static MemoryMappedFile _genomeFile;
        private readonly string _twoBitFile;
        private readonly uint _reserved;
        public uint Reserved { get { return _reserved; }}
        private string FileError { get { return _twoBitFile ?? _genomeFile.ToString(); } }
        private const uint ForwardSignature = 0x1A412743;
        private const uint ReverseSignature = 0x4327411A;
        private const uint TwoBitVersion = 0;

        private const int MaxByte = byte.MaxValue + 1;
        private static readonly byte[] EndiannessIndexedByteArray = new byte[MaxByte];
        private static readonly Dictionary<byte, string> ByteToTetraNucleotideMap = new Dictionary<byte, string>(MaxByte); 

        public const uint NucsPerByte = 4;
        
        internal static string ByteToTetraNucleotide(byte value)
        {
            return ByteToTetraNucleotideMap[EndiannessIndexedByteArray[value]];
        }

        public TwoBitGenomeReader(FileInfo file) : this(file.FullName) { }
        public TwoBitGenomeReader(string file)
        {
            _twoBitFile = file;
            int numSequences;
            Tuple<string, uint>[] contigInfos;
            _genomeFile = MemoryMappedFile.CreateFromFile(file, FileMode.Open);
            var bytePosition = 4U;
            using (var br = _genomeFile.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read))
            {
                ReadSignature(br);
                ReadVersion(br, bytePosition);
                bytePosition += 4;
                unchecked { numSequences = (int)ReadUint(br, bytePosition); }
                bytePosition += 4;
                _contigs = new IGenomeContig[numSequences];
                _reserved = ReadUint(br, bytePosition);
                bytePosition += 4;
                if (_reserved != 0)
                    throw new InvalidTwoBitFileException(FileError, "Reserved is something other than 0: " + _reserved);
                contigInfos = new Tuple<string, uint>[numSequences];
                for (var i = 0; i < numSequences; i++)
                    contigInfos[i] = ReadSequenceIndex(br, bytePosition, out bytePosition);
            }
            Parallel.For(0, numSequences, i =>
            {
                _contigs[i] = new TwoBitGenomeContig(contigInfos[i].Item1, contigInfos[i].Item2,
                    i == numSequences - 1
                        ? 0 // means end of file.
                        : (uint)Math.Abs((long)contigInfos[i + 1].Item2 - (contigInfos[i].Item2)));
            });
            Array.Reverse(_contigs); // turns out the order is always reverse.
        }

        private bool ReadSignature(MemoryMappedViewAccessor mmva)
        {
            var signature = ReadUintLiteral(mmva); // only time it will need to be read literally.
            bool reverse = false;
            if (signature == ForwardSignature) reverse = true;
            else if (signature != ReverseSignature)
                throw new InvalidTwoBitFileException(FileError, "Invalid Signature found: " + signature);
            byte i = 0;
            var temparray = new[] { 'T', 'C', 'A', 'G' };
            foreach (var c in temparray)
                foreach (var d in temparray)
                    foreach (var e in temparray)
                        foreach (var f in temparray)
                        {
                            EndiannessIndexedByteArray[i] = reverse ? ReverseByte(i) : i;
                            ByteToTetraNucleotideMap[i++] = new string(new[] { c, d, e, f });
                        }
            return reverse;
        }

        private void ReadVersion(MemoryMappedViewAccessor mmva, uint bytePosition)
        {
            var version = ReadUint(mmva, bytePosition);
            if (version != TwoBitVersion)
                throw new InvalidTwoBitFileException(FileError, "Version is something other than 0: " + version);
        }

        private Tuple<string, uint> ReadSequenceIndex(MemoryMappedViewAccessor mmva, uint start, out uint bytePosition)
        {
            var nameSize = ReadNBytes(mmva, 1, start)[0];
            bytePosition = start + 1;
            var name = new char[nameSize];
            for (var i = 0; i < nameSize; i++)
            {
                name[i] = (char)ReadNBytes(mmva, 1, bytePosition)[0];
                bytePosition++;
            }
            start = bytePosition;
            bytePosition += 4;
            return Tuple.Create(new string(name), ReadUint(mmva, start));
        }

        private static byte[] ReadNBytes(MemoryMappedViewAccessor mmva, uint n, uint start = 0)
        {
            var byteArray = new byte[n];
            for (var i = 0; i < n; i++)
                byteArray[i] = EndiannessIndexedByteArray[mmva.ReadByte(start + i)];
            if (byteArray.Length > 1 && BitConverter.IsLittleEndian) Array.Reverse(byteArray);
            return byteArray;
        }

        private uint ReadUintLiteral(MemoryMappedViewAccessor mmva) { return mmva.ReadUInt32(0); }

        internal static uint ReadUint(MemoryMappedViewAccessor mmva, uint start = 0)
        {
            return BitConverter.ToUInt32(ReadNBytes(mmva, 4, start), 0);
        }

        private byte ReverseByte(byte b)
        {
            int rev = (b >> 4) | ((b & 0xf) << 4);
            rev = ((rev & 0xcc) >> 2) | ((rev & 0x33) << 2);
            rev = ((rev & 0xaa) >> 1) | ((rev & 0x55) << 1);
            return (byte)rev;
        }

        public void Dispose()
        {
            foreach (var contig in _contigs)
                ((TwoBitGenomeContig)contig).Dispose();
            _genomeFile.Dispose();
        }
    }
}
