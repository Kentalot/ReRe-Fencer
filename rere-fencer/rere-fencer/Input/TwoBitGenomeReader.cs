using System;
using System.CodeDom;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.IO;
using System.IO.MemoryMappedFiles;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntervalTreeLib;
using rere_fencer.Exceptions;

namespace rere_fencer.Input
{

    public class TwoBitGenomeReader : IGenomeReader, IDisposable
    {
        private interface ITwoBitGenomeSubSequence : IGenomeContig
        {
            ulong Start { get; }
            ulong End { get; }
        }

        private abstract class TwoBitGenomeSubSequence : ITwoBitGenomeSubSequence
        {
            public string Name { get; protected set; }
            public ulong Start { get; protected set; }
            public ulong End { get; protected set; }
            public ulong Length { get; protected set; }
            public abstract bool ContainsNs { get; }
            public abstract bool ContainsMaskedSequences { get; }

            public TwoBitGenomeSubSequence(string name, ulong start, ulong end)
            { Name = name; Start = start; End = end; Length = end - start + 1; }

            protected abstract string GetSubSequence(ulong start, ulong end, ulong length, 
                bool ignoreMasks = false, bool ignoreNs = false);

            public string GetSequence(ulong start, ulong end, bool ignoreMasks = false, bool ignoreNs = false)
            {
                if (start < Start || end > End) throw new OutOfContigRangeException(Name, start, end, 
                    string.Format("Inside {0}'s GetSequence method.", GetType().Name));
                var length = end - start + 1;
                return GetSubSequence(Start - start + 1, End - end + 1, length, ignoreMasks, ignoreNs);
            }
            
            public char GetNucleotideAt(ulong position)
            {
                return GetSequence(position, position)[0];
            }

            public IEnumerable<char> GetNucleotides(ulong start, ulong end, bool ignoreMasks = false, bool ignoreNs = false)
            {
                return GetSequence(start, end, ignoreMasks, ignoreNs);
            }
        }

        private class NSubSequence : TwoBitGenomeSubSequence
        {
            public override bool ContainsNs { get { return true; } }
            public override bool ContainsMaskedSequences { get { return false; } }

            public NSubSequence(string name, uint start, uint end) : base(name, start, end) { }
            
            protected override string GetSubSequence(ulong start, ulong end, ulong length, 
                bool ignoreMasks = false, bool ignoreNs = false)
            {
                return ignoreNs ? "" : new string('N', (int) length);
            }
        }

        private class NormalSubSequence : TwoBitGenomeSubSequence, IDisposable
        {
            public override bool ContainsNs { get { return false; } }
            public override bool ContainsMaskedSequences { get { return false; } }

            protected ushort LeftNucsToTrim { get; private set; } // these are partial bytes on the left that are
            protected ushort RightNucsToTrim { get; private set; } // from the previous sequence and need to be trimmed
            protected readonly MemoryMappedViewAccessor _sequenceAccessor;

            public NormalSubSequence(string name, uint start, uint end, ushort leftNucsToTrim, ushort rightNucsToTrim, 
                MemoryMappedViewAccessor sequenceAccessor) : base (name, start, end)
            {
                LeftNucsToTrim = leftNucsToTrim;
                RightNucsToTrim = rightNucsToTrim;
                _sequenceAccessor = sequenceAccessor;
            }

            /*public string GetWholeSequence()
            {
                var halfLength = Length/2;
                var leftoverNucs = Length - (halfLength)*2;
                return string.Format("{0}{1}", GetSubSequence(halfLength + leftoverNucs, true), 
                    GetSubSequence(halfLength, false));                
            }*/

            protected override string GetSubSequence(ulong start, ulong end, ulong length, 
                bool ignoreMasks = false, bool ignoreNs = false)//uint length, bool fromLeft)
            {
                var sb = new StringBuilder((int) length);
                var paddedStart = start + LeftNucsToTrim - 1;
                var leftSubStringIndex = (int) paddedStart%4; // calculate the 
                var lastByte = (end + LeftNucsToTrim - 1) / 4;
                /*var nucsLeftover = (int) (paddedLength - bytesToRead*4); // these are partial bytes that need to be added
                if (fromLeft)
                {
                    bytePosition = 0;
                    lastByte = bytesToRead;
                    if (LeftNucsToTrim > 0) //prepend the first nucs and read one less byte.
                        sb.Append(ByteToTetraNucleotide(_sequenceAccessor.ReadByte(bytePosition++))
                                .Substring(LeftNucsToTrim - 1));
                }
                else
                {
                    lastByte = _sequenceAccessor.Capacity;
                    bytePosition = lastByte - bytesToRead;
                    if (nucsLeftover > 0) //prepend the extra nucs.
                        sb.Append(ByteToTetraNucleotide(_sequenceAccessor.ReadByte(bytePosition - 1))
                            .Substring(nucsLeftover - 1));
                    if (RightNucsToTrim > 0) lastByte--; // read one less byte so I can postpend the last nucs
                }*/

                for (var bytePosition = paddedStart / 4; bytePosition <= lastByte; bytePosition++)
                    sb.Append(ByteToTetraNucleotide(_sequenceAccessor.ReadByte(bytePosition)));

                /*if (fromLeft)
                {
                    if (nucsLeftover > 0)
                        sb.Append(ByteToTetraNucleotide(_sequenceAccessor.ReadByte(lastByte)).Substring(0, 4 - nucsLeftover));
                }
                else if (RightNucsToTrim > 0)
                    sb.Append(ByteToTetraNucleotide(_sequenceAccessor.ReadByte(lastByte)).Substring(0, 4 - RightNucsToTrim));*/
                return sb.ToString().Substring(leftSubStringIndex, (int) length);                
            }

            public void Dispose()
            {
                _sequenceAccessor.Dispose();
            }
        }

        private class MaskedSubSequence : NormalSubSequence
        {
            public override bool ContainsMaskedSequences { get { return true; } }

            public MaskedSubSequence(string name, uint start, uint end, ushort leftNucsToTrim, ushort rightNucsToTrim, 
                MemoryMappedViewAccessor sequenceAccessor)
                : base(name, start, end, leftNucsToTrim, rightNucsToTrim, sequenceAccessor) { }

            protected override string GetSubSequence(ulong start, ulong end, ulong length, 
                bool ignoreMasks = false, bool ignoreNs = false)
            {
                var tempstr = base.GetSubSequence(start, end, length, ignoreMasks, ignoreNs);
 	            return ignoreMasks ? tempstr : tempstr.ToLower();
            }
        }

        private class TwoBitGenomeContig : TwoBitGenomeSubSequence
        {
            public override bool ContainsNs { get { return _containsNs; } } 
            private readonly bool _containsNs;
            public override bool ContainsMaskedSequences { get { return _containsMaskedSequences; } } 
            private readonly bool _containsMaskedSequences;

            private readonly List<ITwoBitGenomeSubSequence> _subSequences;

            public TwoBitGenomeContig(string name, ulong length, List<ITwoBitGenomeSubSequence> subSequences)
                : base(name, 1, length)
            {
                _subSequences = subSequences;
                _containsNs = _subSequences.Any(s => s.ContainsNs);
                _containsMaskedSequences = _subSequences.Any(s => s.ContainsMaskedSequences);
            }

            protected override string GetSubSequence(ulong start, ulong end, ulong length, 
                bool ignoreMasks = false, bool ignoreNs = false)
            {
                var i = BinarySearchOverlappingRegions(start, 0, _subSequences.Count - 1);
                if (end < _subSequences[i].End) // all contained in one subSequence;
                    return _subSequences[i].GetSequence(start, end, ignoreMasks, ignoreNs);
                var endingIndex = BinarySearchOverlappingRegions(end, i + 1, _subSequences.Count - 1);
                var seqString = _subSequences[i].GetSequence(start, _subSequences[i++].End, ignoreMasks, ignoreNs);
                for (; i < endingIndex; i++)
                    seqString += _subSequences[i].GetSequence(_subSequences[i].Start, _subSequences[i].End, ignoreMasks, ignoreNs);
                seqString += _subSequences[i].GetSequence(_subSequences[i].Start, end, ignoreMasks, ignoreNs);
                return seqString;
            }

            private int BinarySearchOverlappingRegions(ulong position, int begin, int end)
            {
                if (begin == end) return begin; // already found it.
                var middle = (begin + end)/2;
                if (position < _subSequences[middle].Start) // before middle
                    return BinarySearchOverlappingRegions(position, begin, middle - 1);
                if (position > _subSequences[middle].End) // after middle
                    return BinarySearchOverlappingRegions(position, middle + 1, end);
                return middle;
            }
        }

        public bool IsZeroBasedCoordinates { get { return false; } }
        public bool SupportsMasking { get { return true; } }
        public bool SupportsNs { get { return true; } }
        public bool SupportsIupacAmbiguityCodes { get { return false; } }

        private readonly List<IGenomeContig> _contigs;
        public List<IGenomeContig> Contigs { get { return _contigs; } }

        private readonly MemoryMappedFile _mmf;
        private readonly ITwoBitGenomeSubSequence[] _sequenceInfos;
        private readonly bool _reverse = false;
        private string _twoBitFile = null;
        private readonly uint _reserved;
        public uint Reserved { get { return _reserved; }}
        private string _fileError { get { return _twoBitFile ?? _mmf.ToString(); } }
        private const uint ForwardSignature = 0x1A412743;
        private const uint ReverseSignature = 0x4327411A;
        private const uint TwoBitVersion = 0;

        internal static readonly byte[] EndiannessIndexedByteArray = new byte[256];

        private class SequenceInfo
        {
            private readonly string _sequenceName;
            public string SequenceName { get { return _sequenceName; } }
            private readonly uint _origin;
            public uint Origin { get { return _origin; } }
            private readonly SequenceData _data;

            public SequenceInfo(string sequenceName, uint origin, BinaryReader br)
            {
                _sequenceName = sequenceName;
                _origin = origin;
            }

            public void LoadSequenceData { }
        }

        //private static char[] bit_chars = {'T', 'C', 'A', 'G'};

        public enum TetraNucleotide : byte
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
        }
        
        internal static string ByteToTetraNucleotide(byte value)
        {
            return ((TetraNucleotide) EndiannessIndexedByteArray[value]).ToString();
        }

        public TwoBitGenomeReader(string file)
        {
            _twoBitFile = file;
            uint numSequences;
            using (var stream = File.OpenRead(file))
            {
                using (var br = new BinaryReader(stream))
                {
                    ReadSignature(br);
                    ReadVersion(br);
                    numSequences = ReadUint(br);
                    //_sequenceInfos = new SequenceInfo[ReadUint(br)];
                    _reserved = ReadUint(br);
                }
            }
            _mmf = MemoryMappedFile.CreateFromFile(file, FileMode.Open);
            ReadSequenceIndex();
            
        }

        private void ReadSequenceIndex()
        {
            throw new NotImplementedException();
        }

        public TwoBitGenomeReader(FileInfo file) : this(file.FullName) { }

        private bool ReadSignature(BinaryReader br)
        {
            var signature = ReadUint(br);
            bool reverse = false;
            if (signature == ForwardSignature) reverse = true;
            else if (signature != ReverseSignature) 
                throw new InvalidTwoBitFileException(_fileError, "Invalid Signature found: " + signature);
            for (byte i = 0; i < EndiannessIndexedByteArray.Length; i++)
                EndiannessIndexedByteArray[i] = reverse ? ReverseByte(i) : i;
            return reverse;
        }

        private void ReadVersion(BinaryReader br)
        {
            var version = ReadUint(br);
            if (version != TwoBitVersion) throw new InvalidTwoBitFileException(_fileError, "Version is something other than 0: " + version);
        }

        private byte[] ReadNBytes(BinaryReader br, int n)
        {
            var byteArray = new byte[n];
            for (var i = 0; i < n; i++)
                byteArray[i] = EndiannessIndexedByteArray[br.ReadByte()];
            return byteArray;
        }

        private uint ReadUint(BinaryReader br, bool modify = false)
        {
            return modify ? BitConverter.ToUInt32(ReadNBytes(br, 4), 0) : br.ReadUInt32();
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
            _mmf.Dispose();
        }

    }
}
