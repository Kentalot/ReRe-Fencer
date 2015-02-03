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
        private class TwoBitGenomeContig : IGenomeContig
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
            public IEnumerable<IContigRegion> NormalRegions
            { 
                get { return null; } 
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

            public IEnumerable<char> GetSequence(uint start, uint end,
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                if (start < 1 || end > Length) throw new OutOfContigRangeException(Name, start, end,
                    string.Format("Inside {0}'s GetSequence method.", GetType().Name));
                if (end < start) throw new ArgumentException("Sequence end was less than start, what's up with that?");
                var length = end - start + 1;
                uint charPosition = 0, bytePosition = start/NucsPerByte;
                var leftStartIndex = start % NucsPerByte;
                var numBytes = (int) Math.Ceiling((length + leftStartIndex)/(decimal)NucsPerByte);
                var overlappingNBlocks = skipNs ? new List<int>() : GetAllOverlappingBlocks(start, end, true);
                var overlappingMaskBlocks = ignoreMasks ? new List<int>() : GetAllOverlappingBlocks(start, end, false);
                if (!overlappingNBlocks.Any()) skipNs = true;
                if (!overlappingMaskBlocks.Any()) ignoreMasks = true;
                int nIndex = 0, maskIndex = 0;
                var position = start;
                var inNBlock = !skipNs &&
                               IsOverlapping(position, overlappingNBlocks, _nBlockStarts, _nBlockEnds, nIndex);
                var inMaskBlock = !ignoreMasks &&
                                  IsOverlapping(position, overlappingMaskBlocks, _maskBlockStarts, _maskBlockEnds, maskIndex);
                for (var j = 0; j < numBytes; j++)
                {
                    var tetNuc = ByteToTetraNucleotide(_contigAccessor.ReadByte(bytePosition));
                    for (var i = j == 0 ? (int) leftStartIndex : 0; i < NucsPerByte; i++)
                    {
                        if (inNBlock)
                        {
                            yield return 'N';
                            charPosition++;
                            position++;
                            if (charPosition == length) yield break;
                            inNBlock = IsOverlapping(position, overlappingNBlocks, _nBlockStarts, _nBlockEnds, nIndex);
                            if (inNBlock) continue;
                            if (++nIndex == overlappingNBlocks.Count)
                                skipNs = true; // already at the end of the nBlocks.
                            else
                                inNBlock = IsOverlapping(position, overlappingNBlocks, _nBlockStarts, _nBlockEnds, nIndex);
                        }
                        else if (inMaskBlock)
                        {
                            yield return char.ToLower(tetNuc[i], CultureInfo.CurrentCulture);
                            charPosition++;
                            position++;
                            if (charPosition == length) yield break;
                            inMaskBlock = IsOverlapping(position, overlappingMaskBlocks, _maskBlockStarts, _maskBlockEnds, maskIndex);
                            if (inMaskBlock) continue;
                            if (++maskIndex == overlappingMaskBlocks.Count) ignoreMasks = true;
                            else
                                inMaskBlock =
                                  IsOverlapping(position, overlappingMaskBlocks, _maskBlockStarts, _maskBlockEnds, maskIndex);
                        }
                        else
                        {
                            yield return tetNuc[i];
                            charPosition++;
                            position++;
                            if (charPosition == length) yield break;
                            if (!ignoreMasks)
                                inMaskBlock =
                                  IsOverlapping(position, overlappingMaskBlocks, _maskBlockStarts, _maskBlockEnds, maskIndex);
                            if (!skipNs && !inMaskBlock)
                                inNBlock = IsOverlapping(position, overlappingNBlocks, _nBlockStarts, _nBlockEnds, nIndex);
                        }
                    }
                    bytePosition++;
                }
            }

            private bool IsOverlapping(uint position, IList<int> overlappingBlocks, uint[] blockStarts, uint[] blockEnds, int index)
            {
                return position >= blockStarts[overlappingBlocks[index]] && position <= blockEnds[overlappingBlocks[index]];
            }

            private IList<int> GetAllOverlappingBlocks(uint start, uint end, bool searchNs)
            {
                var starts = searchNs ? _nBlockStarts : _maskBlockStarts;
                var ends = searchNs ? _nBlockEnds : _maskBlockEnds;
                var tempList = new List<int>();
                if (!starts.Any()) return tempList;
                var firstPossibleOverlap = BinarySearchForFirstOverlappingBlock(starts, ends, start, 0, starts.Length);
                if (firstPossibleOverlap < 0 || end < ends[firstPossibleOverlap]) // no overlaps possible.
                    return tempList;
                var lastPossibleOverlap = BinarySearchForLastOverlappingBlock(starts, ends, end, firstPossibleOverlap, starts.Length);
                if (lastPossibleOverlap < 0) return new List<int>{firstPossibleOverlap};
                for (var i = firstPossibleOverlap; i <= lastPossibleOverlap; i++)
                    tempList.Add(i);
                return tempList;
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

            public IEnumerable<char> GetNucleotides(uint start, uint end, bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                return GetSequence(start, end, ignoreMasks, skipMasks, skipNs);
            }

            private uint GetBlockStartsAndEnds(MemoryMappedViewAccessor tempViewer, uint start, uint[] blockStarts, uint[] blockEnds)
            {
                const uint uintByteCount = 4;
                var blockCount = blockStarts.Length;
                 // not sure if this will break if uint > Int32.MaxValue
                for (var i = 0; i < blockCount; i++)
                {
                    blockStarts[i] = ReadUint(tempViewer, start);
                    start += uintByteCount;
                }
                for (var i = 0; i < blockCount; i++)
                {
                    blockEnds[i] = blockStarts[i] + ReadUint(tempViewer, start) - 1;
                    start += uintByteCount;
                }
                return start;
            }
        }

        public bool IsZeroBasedCoordinates { get { return true; } }
        public bool SupportsMasking { get { return true; } }
        public bool SupportsNs { get { return true; } }
        public bool SupportsIupacAmbiguityCodes { get { return false; } }

        private readonly IGenomeContig[] _contigs;
        public IReadOnlyList<IGenomeContig> Contigs { get { return Array.AsReadOnly(_contigs); } }

        private static MemoryMappedFile _genomeFile;
        private static MemoryMappedViewAccessor[] _genomeAccessors;
        private readonly bool _reverse = false;
        private string _twoBitFile;
        private readonly uint _reserved;
        public uint Reserved { get { return _reserved; }}
        private string _fileError { get { return _twoBitFile ?? _genomeFile.ToString(); } }
        private const uint ForwardSignature = 0x1A412743;
        private const uint ReverseSignature = 0x4327411A;
        private const uint TwoBitVersion = 0;

        private const int maxByte = byte.MaxValue + 1;
        private static readonly byte[] EndiannessIndexedByteArray = new byte[maxByte];
        private static readonly Dictionary<byte, string> ByteToTetraNucleotideMap = new Dictionary<byte, string>(maxByte); 
        private static readonly Dictionary<byte, string> ByteToMaskedTetraNucleotideMap = new Dictionary<byte, string>(maxByte); 
        
        //private static char[] bit_chars = {'T', 'C', 'A', 'G'};

        public const uint NucsPerByte = 4;
        /*public enum TetraNucleotide : byte
        {
            TTTT = 0x00, TTTC, TTTA, TTTG, TTCT, TTCC, TTCA, TTCG, TTAT, TTAC, TTAA, TTAG, TTGT, TTGC, TTGA, TTGG,
            TCTT, TCTC, TCTA, TCTG, TCCT, TCCC, TCCA, TCCG, TCAT, TCAC, TCAA, TCAG, TCGT, TCGC, TCGA, TCGG,
            TATT, TATC, TATA, TATG, TACT, TACC, TACA, TACG, TAAT, TAAC, TAAA, TAAG, TAGT, TAGC, TAGA, TAGG,
            TGTT, TGTC, TGTA, TGTG, TGCT, TGCC, TGCA, TGCG, TGAT, TGAC, TGAA, TGAG, TGGT, TGGC, TGGA, TGGG,
            CTTT, CTTC, CTTA, CTTG, CTCT, CTCC, CTCA, CTCG, CTAT, CTAC, CTAA, CTAG, CTGT, CTGC, CTGA, CTGG,
            CCTT, CCTC, CCTA, CCTG, CCCT, CCCC, CCCA, CCCG, CCAT, CCAC, CCAA, CCAG, CCGT, CCGC, CCGA, CCGG,
            CATT, CATC, CATA, CATG, CACT, CACC, CACA, CACG, CAAT, CAAC, CAAA, CAAG, CAGT, CAGC, CAGA, CAGG,
            CGTT, CGTC, CGTA, CGTG, CGCT, CGCC, CGCA, CGCG, CGAT, CGAC, CGAA, CGAG, CGGT, CGGC, CGGA, CGGG,
            ATTT, ATTC, ATTA, ATTG, ATCT, ATCC, ATCA, ATCG, ATAT, ATAC, ATAA, ATAG, ATGT, ATGC, ATGA, ATGG,
            ACTT, ACTC, ACTA, ACTG, ACCT, ACCC, ACCA, ACCG, ACAT, ACAC, ACAA, ACAG, ACGT, ACGC, ACGA, ACGG,
            AATT, AATC, AATA, AATG, AACT, AACC, AACA, AACG, AAAT, AAAC, AAAA, AAAG, AAGT, AAGC, AAGA, AAGG,
            AGTT, AGTC, AGTA, AGTG, AGCT, AGCC, AGCA, AGCG, AGAT, AGAC, AGAA, AGAG, AGGT, AGGC, AGGA, AGGG,
            GTTT, GTTC, GTTA, GTTG, GTCT, GTCC, GTCA, GTCG, GTAT, GTAC, GTAA, GTAG, GTGT, GTGC, GTGA, GTGG,
            GCTT, GCTC, GCTA, GCTG, GCCT, GCCC, GCCA, GCCG, GCAT, GCAC, GCAA, GCAG, GCGT, GCGC, GCGA, GCGG,
            GATT, GATC, GATA, GATG, GACT, GACC, GACA, GACG, GAAT, GAAC, GAAA, GAAG, GAGT, GAGC, GAGA, GAGG,
            GGTT, GGTC, GGTA, GGTG, GGCT, GGCC, GGCA, GGCG, GGAT, GGAC, GGAA, GGAG, GGGT, GGGC, GGGA, GGGG
        }*/
        
        internal static string ByteToTetraNucleotide(byte value, bool masked = false)
        {
            return masked
                ? ByteToMaskedTetraNucleotideMap[EndiannessIndexedByteArray[value]]
                : ByteToTetraNucleotideMap[EndiannessIndexedByteArray[value]];
        }

        public TwoBitGenomeReader(FileInfo file) : this(file.FullName) { }
        public TwoBitGenomeReader(string file)
        {
            _twoBitFile = file;
            int numSequences;
            Tuple<string, uint>[] contigInfos;
            using (var stream = File.OpenRead(file))
            {
                using (var br = new BinaryReader(stream))
                {
                    ReadSignature(br);
                    ReadVersion(br);
                    unchecked { numSequences = (int) ReadUint(br); }
                    _contigs = new IGenomeContig[numSequences];
                    _genomeAccessors = new MemoryMappedViewAccessor[numSequences];
                    _reserved = ReadUint(br);
                    if (_reserved != 0)
                        throw new InvalidTwoBitFileException(_fileError, "Reserved is something other than 0: " + _reserved);
                    contigInfos = new Tuple<string, uint>[numSequences];
                    for (var i = 0; i < numSequences; i++)
                        contigInfos[i] = ReadSequenceIndex(br);
                }
            }
            _genomeFile = MemoryMappedFile.CreateFromFile(file, FileMode.Open);
            for (var i = 0; i < numSequences; i++)
                _contigs[i] = new TwoBitGenomeContig(contigInfos[i].Item1, contigInfos[i].Item2,
                    i == numSequences - 1 
                        ? 0 // means end of file.
                        : (uint) Math.Abs((long) contigInfos[i + 1].Item2 - (contigInfos[i].Item2)));
                    //GetSequenceContent(contigInfos[i].Item1, contigInfos[i].Item2, i == numSequences - 1 ? 0 : contigInfos[i + 1].Item2, i);
            Array.Reverse(_contigs); // turns out the order is always reverse.
        }

        private bool ReadSignature(BinaryReader br)
        {
            var signature = ReadUint(br, true); // only time it will need to be read literally.
            bool reverse = false;
            if (signature == ForwardSignature) reverse = true;
            else if (signature != ReverseSignature) 
                throw new InvalidTwoBitFileException(_fileError, "Invalid Signature found: " + signature);
            //for (byte i = 0; i < EndiannessIndexedByteArray.Length; i++)
            //    EndiannessIndexedByteArray[i] = reverse ? ReverseByte(i) : i;
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
            temparray = new[] { 't', 'c', 'a', 'g' };
            i = 0;
            foreach (var c in temparray)
                foreach (var d in temparray)
                    foreach (var e in temparray)
                        foreach (var f in temparray)
                            ByteToMaskedTetraNucleotideMap[i++] = new string(new[] { c, d, e, f });
            return reverse;
        }

        private void ReadVersion(BinaryReader br)
        {
            var version = ReadUint(br);
            if (version != TwoBitVersion) 
                throw new InvalidTwoBitFileException(_fileError, "Version is something other than 0: " + version);
        }

        private Tuple<string, uint> ReadSequenceIndex(BinaryReader br)
        {
            var nameSize = ReadNBytes(br, 1)[0];
            var name = new char[nameSize];
            for (var i = 0; i < nameSize; i++)
                name[i] = (char)ReadNBytes(br, 1)[0];
            return Tuple.Create(new string(name), ReadUint(br));
        }

        private byte[] ReadNBytes(BinaryReader br, uint n)
        {
            var byteArray = new byte[n];
            for (var i = 0; i < n; i ++)
                byteArray[i] = EndiannessIndexedByteArray[br.ReadByte()];
            if (BitConverter.IsLittleEndian) Array.Reverse(byteArray);
            return byteArray;
        }

        private static byte[] ReadNBytes(MemoryMappedViewAccessor mmva, uint n, uint start = 0)
        {
            var byteArray = new byte[n];
            for (var i = 0; i < n; i++)
                byteArray[i] = EndiannessIndexedByteArray[mmva.ReadByte(start + i)];
            if (BitConverter.IsLittleEndian) Array.Reverse(byteArray);
            return byteArray;
        }

        private uint ReadUint(BinaryReader br, bool readLiteral = false)
        {
            return readLiteral ? br.ReadUInt32() : BitConverter.ToUInt32(ReadNBytes(br, 4), 0);
        }

        internal static uint ReadUint(MemoryMappedViewAccessor mmva, uint start = 0)
        {
            return BitConverter.ToUInt32(ReadNBytes(mmva, 4, start), 0);
        }

        internal static byte ReverseByte(byte b)
        {
            int rev = (b >> 4) | ((b & 0xf) << 4);
            rev = ((rev & 0xcc) >> 2) | ((rev & 0x33) << 2);
            rev = ((rev & 0xaa) >> 1) | ((rev & 0x55) << 1);
            return (byte)rev; 
        }

        public void Dispose()
        {
            foreach (var mmva in _genomeAccessors)
                if (mmva != null) mmva.Dispose();
            _genomeFile.Dispose();
        }

    }
}
