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
        private interface ITwoBitGenomeSubcontig
        {
            string Name { get; }
            IContigRegion RegionInfo { get; }

            IEnumerable<char> GetSequence(uint start, uint end, bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false);
        }

        private abstract class TwoBitGenomeSubcontig : ITwoBitGenomeSubcontig
        {
            public string Name { get; private set; }
            public IContigRegion RegionInfo { get; private set; }

            public TwoBitGenomeSubcontig(string name, IContigRegion region)
            { 
                Name = name; RegionInfo = region;
            }

            protected abstract IEnumerable<char> GetSubSequence(uint start, uint end, uint length, 
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false);

            public IEnumerable<char> GetSequence(uint start, uint end,
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                if (start < RegionInfo.Start || end > RegionInfo.End) throw new OutOfContigRangeException(Name, start, end, 
                    string.Format("Inside {0}'s GetSequence method.", GetType().Name));
                if (end < start) throw new ArgumentException("Sequence end was less than start, what's up with that?");
                var length = end - start + 1;
                return GetSubSequence(start, end, length, ignoreMasks, skipMasks, skipNs);
            }
        }

        private class NormalSubcontig : TwoBitGenomeSubcontig
        {
            private readonly ushort _leftNucsToTrim; // these are partial bytes on the left that are trimmed before returning
            private readonly MemoryMappedViewAccessor _sequenceAccessor;

            public NormalSubcontig(string name, IContigRegion region, ushort leftNucsToTrim, 
                MemoryMappedViewAccessor sequenceAccessor) 
                : base (name, region)
            {
                _leftNucsToTrim = leftNucsToTrim;
                _sequenceAccessor = sequenceAccessor;
            }
            
            protected override IEnumerable<char> GetSubSequence(uint start, uint end, uint length, 
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)//uint length, bool fromLeft)
            {
                var intLength = (int) length;
                var returnChars = new char[intLength];
                var bytePosition = start / NucsPerByte; // calculate the starting position to read from.
                var leftSubStringIndex = (int)(start % NucsPerByte); // calculate the leftovers
                //var lastByte = bytePosition + (uint) Math.Ceiling(paddedLength/(double)NucsPerByte);
                var lastByte = (int) Math.Ceiling((end + (double)_leftNucsToTrim)/ NucsPerByte);

                var charPosition = 0;
                Console.WriteLine("Reading from {0} to {1} starting from BytePostion {2} and lastByte is {3}", start,
                    end, bytePosition, lastByte);
                var tetNuc = ByteToTetraNucleotide(_sequenceAccessor.ReadByte(bytePosition++));
                for (var i = leftSubStringIndex; i < NucsPerByte && charPosition < returnChars.Length; i++)
                {
                    yield return tetNuc[i];
                    charPosition++;
                }

                for (; bytePosition <= lastByte; bytePosition++)
                {
                    tetNuc = ByteToTetraNucleotide(_sequenceAccessor.ReadByte(bytePosition));
                    for (var i = 0; i < NucsPerByte && charPosition < returnChars.Length; i++)
                        // charPos < returnChars.Length short circuits
                    {
                        yield return tetNuc[i];
                        charPosition++;
                    }
                }
            }
        }

        private class MaskedSubcontig : NormalSubcontig
        {
            public MaskedSubcontig(string name, IContigRegion region, ushort leftNucsToTrim, //ushort rightNucsToTrim, 
                MemoryMappedViewAccessor sequenceAccessor)
                : base(name, region, leftNucsToTrim, sequenceAccessor) { }

            protected override IEnumerable<char> GetSubSequence(uint start, uint end, uint length, 
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                if (skipMasks) return String.Empty;
                var tempstr = base.GetSubSequence(start, end, length, ignoreMasks, skipMasks, skipNs);
 	            return ignoreMasks ? tempstr : tempstr.Select(c => char.ToLower(c, CultureInfo.CurrentCulture));
            }
        }

        private class NSubcontig : TwoBitGenomeSubcontig
        {
            public NSubcontig(string name, IContigRegion region) : base(name, region) { }

            protected override IEnumerable<char> GetSubSequence(uint start, uint end, uint length,
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                return skipNs ? "" : new string('N', (int)length);
            }
        }

        private class TwoBitGenomeContig : IGenomeContig
        {
            public string Name { get; private set; }
            public uint Length { get; private set; }
            public bool ContainsNs { get { return NRegions.Any(); } }
            public bool ContainsMaskedSequences { get { return MaskedRegions.Any(); } }
            public IEnumerable<IContigRegion> NRegions 
                { get { return _subcontigs.Where(s => s.GetType() == typeof(NSubcontig)).Select(s => s.RegionInfo); } }
            public IEnumerable<IContigRegion> MaskedRegions
                { get { return _subcontigs.Where(s => s.GetType() == typeof(MaskedSubcontig)).Select(s => s.RegionInfo); } }
            public IEnumerable<IContigRegion> NormalRegions
            { 
                get
                {
                    return _subcontigs.Where(s =>
                    {
                        var type = s.GetType();
                        return type != typeof (NSubcontig) && type != typeof (MaskedSubcontig);
                    }).Select(s => s.RegionInfo);
                } 
            }

            private readonly List<ITwoBitGenomeSubcontig> _subcontigs;

            public TwoBitGenomeContig(string name, uint length, List<ITwoBitGenomeSubcontig> subcontigs)
            {
                Name = name;
                Length = length;
                _subcontigs = subcontigs;
            }

            public IEnumerable<char> GetSequence(uint start, uint end,
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                if (start < 1 || end > Length) throw new OutOfContigRangeException(Name, start, end,
                    string.Format("Inside {0}'s GetSequence method.", GetType().Name));
                if (end < start) throw new ArgumentException("Sequence end was less than start, what's up with that?");
                var length = end - start + 1;
                var charPosition = 0;
                int position, endingIndex;

                if (Length - end > start - 1) // end is closer to the middle, so do the first search using end
                {
                    endingIndex = BinarySearchForIndex(end, 0, _subcontigs.Count - 1);
                    if (start >= _subcontigs[endingIndex].RegionInfo.Start) // all contained in one subSequence;
                    {
                        foreach (var returnChar in _subcontigs[endingIndex].GetSequence(start, end, ignoreMasks, skipMasks, skipNs))
                        {
                            yield return returnChar;
                            charPosition++;
                            if (charPosition > length) yield break;
                        }
                    }
                    position = BinarySearchForIndex(start, 0, endingIndex - 1);
                }
                else
                {
                    position = BinarySearchForIndex(start, 0, _subcontigs.Count - 1);
                    if (end <= _subcontigs[position].RegionInfo.End) // all contained in one subSequence;
                    {
                        foreach (var returnChar in _subcontigs[position].GetSequence(start, end, ignoreMasks, skipMasks, skipNs))
                        {
                            yield return returnChar;
                            charPosition++;
                            if (charPosition > length) yield break;
                        }
                    }
                    endingIndex = BinarySearchForIndex(end, position + 1, _subcontigs.Count - 1);
                }

                foreach (var returnChar in _subcontigs[position].GetSequence(start, _subcontigs[position].RegionInfo.End, ignoreMasks,skipMasks, skipNs))
                {
                    yield return returnChar;
                    charPosition++;
                    if (charPosition > length) yield break;
                }
                position++;

                for (; position < endingIndex; position++)
                {
                    foreach (var returnChar in _subcontigs[position].GetSequence(_subcontigs[position].RegionInfo.Start,
                        _subcontigs[position].RegionInfo.End, ignoreMasks, skipMasks, skipNs))
                    {
                        yield return returnChar;
                        charPosition++;
                        if (charPosition > length) yield break;
                    }
                    
                }

                foreach (var returnChar in _subcontigs[position].GetSequence(_subcontigs[position].RegionInfo.Start, end, ignoreMasks, skipMasks, skipNs))
                {
                    yield return returnChar;
                    charPosition++;
                    if (charPosition > length) yield break;
                }
            }

            private int BinarySearchForIndex(uint position, int begin, int end)
            {
                if (begin == end) return begin; // already found it.
                var middle = (begin + end) / 2;
                if (position < _subcontigs[middle].RegionInfo.Start) // before middle
                    return BinarySearchForIndex(position, begin, middle - 1);
                if (position > _subcontigs[middle].RegionInfo.End) // after middle
                    return BinarySearchForIndex(position, middle + 1, end);
                return middle;
            }

            public char GetNucleotideAt(uint position)
            {
                return GetSequence(position, position).First();
            }

            public IEnumerable<char> GetNucleotides(uint start, uint end, bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                return GetSequence(start, end, ignoreMasks, skipMasks, skipNs);
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
        private readonly uint[] _contigReserved;
        public uint[] ContigReserved { get { return _contigReserved; } }
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
                    _contigReserved = new uint[numSequences];
                    _genomeAccessors = new MemoryMappedViewAccessor[numSequences];
                    _reserved = ReadUint(br);
                    if (_reserved != 0)
                        throw new InvalidTwoBitFileException(_fileError, "Reserved is something other than 0: " + _reserved);
                    contigInfos = new Tuple<string, uint>[numSequences];
                    for (var i = 0; i < numSequences; i++)
                        contigInfos[i] = ReadSequenceIndex(br);
                    Array.Reverse(contigInfos);
                }
            }
            _genomeFile = MemoryMappedFile.CreateFromFile(file, FileMode.Open);
            for (var i = 0; i < numSequences; i++)
                _contigs[i] = GetSequenceContent(contigInfos[i].Item1, contigInfos[i].Item2, i == numSequences - 1 ? 0 : contigInfos[i + 1].Item2, i);
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

        private IGenomeContig GetSequenceContent(string contigName, uint memoryOffset, uint nextContigOffset, int contigIndex)
        {
            var start = 4U;
            SortedList<uint, uint> nBlocks, maskBlocks;
            uint dnaSize;
            using (var tempViewer = _genomeFile.CreateViewAccessor(memoryOffset, 0, MemoryMappedFileAccess.Read))
            {
                dnaSize = ReadUint(tempViewer);
                start = GetBlockCountFromSequenceIndex(tempViewer, start, out nBlocks);
                start = GetBlockCountFromSequenceIndex(tempViewer, start, out maskBlocks);
                _contigReserved[contigIndex] = ReadUint(tempViewer, start);
                start += 4;
            }

            var subcontigs = new List<ITwoBitGenomeSubcontig>(2 * (nBlocks.Count + maskBlocks.Count + 1));
            var byteOffset = start + memoryOffset;
            _genomeAccessors[contigIndex] = _genomeFile.CreateViewAccessor(byteOffset, nextContigOffset, MemoryMappedFileAccess.Read);
            ushort firstIndex = 0;
            uint nextBlock = 0, currentSequencePosition = 0;
            for (int currentNIndex = 0, currentMaskIndex = 0; currentNIndex < nBlocks.Count || currentMaskIndex < maskBlocks.Count; )
            {
                if (currentNIndex < nBlocks.Count && currentSequencePosition == nBlocks.Keys[currentNIndex])
                {
                    nextBlock = currentSequencePosition + nBlocks.Values[currentNIndex++];
                    subcontigs.Add(GenerateNewSubcontig(contigName, currentSequencePosition, nextBlock - 1));
                }
                else if (currentMaskIndex < maskBlocks.Count && currentSequencePosition == maskBlocks.Keys[currentMaskIndex])
                {
                    nextBlock = currentSequencePosition + maskBlocks.Values[currentMaskIndex++];
                    subcontigs.Add(GenerateNewSubcontig(contigName, currentSequencePosition, nextBlock - 1,
                        _genomeAccessors[contigIndex], firstIndex, true));
                }
                else
                {
                    var nBlockKey = currentNIndex == nBlocks.Count ? dnaSize - 1 : nBlocks.Keys[currentNIndex];
                    var maskBlockKey = currentMaskIndex == maskBlocks.Count ? dnaSize - 1 : maskBlocks.Keys[currentMaskIndex];
                    if (currentSequencePosition < nBlockKey && currentSequencePosition < maskBlockKey)
                    {
                        nextBlock = maskBlockKey < nBlockKey ? maskBlockKey : nBlockKey;
                        subcontigs.Add(GenerateNewSubcontig(contigName, currentSequencePosition, nextBlock - 1,
                            _genomeAccessors[contigIndex], firstIndex));
                    }
                }
                currentSequencePosition = nextBlock;
                firstIndex = (ushort)(nextBlock % NucsPerByte);
                //currentBytePosition = nextBlock / NucsPerByte;
            }

            var lastOne = subcontigs.LastOrDefault();
            if (lastOne == null || lastOne.RegionInfo.End != dnaSize) // means that there's one more normal block left to make
                subcontigs.Add(GenerateNewSubcontig(contigName, currentSequencePosition + 1, dnaSize,
                    _genomeAccessors[contigIndex], firstIndex));

            return new TwoBitGenomeContig(contigName, dnaSize, subcontigs);
        }

        private uint GetBlockCountFromSequenceIndex(MemoryMappedViewAccessor tempViewer, uint start, out SortedList<uint, uint> blocks)
        {
            const uint uintByteCount = 4;
            var blockCount = ReadUint(tempViewer, start);
            var position = uintByteCount + start;
            var capacity = blockCount > (uint)Int32.MaxValue ? Int32.MaxValue : (int)blockCount;
            var blockStarts = new List<uint>(capacity); // not sure if this will break if uint > Int32.MaxValue
            for (var i = 0; i < blockCount; i++)
            {
                blockStarts.Add(ReadUint(tempViewer, position));
                position += uintByteCount;
            }
            blocks = new SortedList<uint, uint>(capacity);
            for (var i = 0; i < blockCount; i++)
            {
                blocks.Add(blockStarts[i], ReadUint(tempViewer, position));
                position += uintByteCount;
            }
            return position;
        }

        private ITwoBitGenomeSubcontig GenerateNewSubcontig(string name, uint start, uint end,
            MemoryMappedViewAccessor mmva = null, ushort nucsToLeft = 0, bool isMaskedContig = false)
        {
            var region = new ContigRegion(start, end);
            var newName = string.Format("{0}:{1}-{2}", name, start, end);
            if (mmva == null) return new NSubcontig(newName, region);
            //var numBytes = (uint)Math.Ceiling(end / (decimal)NucsPerByte) - currentBytePosition;
            return isMaskedContig
                ? new MaskedSubcontig(newName, region, nucsToLeft, mmva)
                : new NormalSubcontig(newName, region, nucsToLeft, mmva);
        }

        private byte[] ReadNBytes(BinaryReader br, uint n)
        {
            var byteArray = new byte[n];
            for (var i = 0; i < n; i ++)
                byteArray[i] = EndiannessIndexedByteArray[br.ReadByte()];
            if (BitConverter.IsLittleEndian) Array.Reverse(byteArray);
            return byteArray;
        }

        private byte[] ReadNBytes(MemoryMappedViewAccessor mmva, uint n, uint start = 0)
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

        private uint ReadUint(MemoryMappedViewAccessor mmva, uint start = 0)
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
                mmva.Dispose();
            _genomeFile.Dispose();
        }

    }
}
