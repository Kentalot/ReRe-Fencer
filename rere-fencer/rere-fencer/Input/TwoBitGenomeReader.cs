using System;
using System.CodeDom;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.IO;
using System.IO.MemoryMappedFiles;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using rere_fencer.Exceptions;

namespace rere_fencer.Input
{

    public class TwoBitGenomeReader : IGenomeReader, IDisposable
    {
        private interface ITwoBitGenomeSubcontig
        {
            string Name { get; }
            uint Start { get; }
            uint End { get; }

            string GetSequence(uint start, uint end, bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false);
        }

        private abstract class TwoBitGenomeSubcontig : ITwoBitGenomeSubcontig
        {
            protected const int NucsPerByte = 4;
            public string Name { get; private set; }
            public uint Start { get; private set; }
            public uint End { get; private set; }

            public TwoBitGenomeSubcontig(string name, uint start, uint end)
            { Name = name; Start = start; End = end; }

            protected abstract string GetSubSequence(uint start, uint end, uint length, 
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false);

            public string GetSequence(uint start, uint end,
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                if (start < Start || end > End) throw new OutOfContigRangeException(Name, start, end, 
                    string.Format("Inside {0}'s GetSequence method.", GetType().Name));
                if (end < start) throw new ArgumentException("Sequence end was less than start, what's up with that?");
                var length = end - start + 1;
                return GetSubSequence(Start - start + 1, End - end + 1, length, ignoreMasks, skipMasks, skipNs);
            }
        }

        private class NormalSubcontig : TwoBitGenomeSubcontig, IDisposable
        {
            private ushort LeftNucsToTrim { get; set; } // these are partial bytes on the left that are
            //private ushort RightNucsToTrim { get; set; } // from the previous sequence and need to be trimmed, seems not necessary
            private readonly MemoryMappedViewAccessor _sequenceAccessor;

            public NormalSubcontig(string name, uint start, uint end, ushort leftNucsToTrim, //ushort rightNucsToTrim, 
                MemoryMappedViewAccessor sequenceAccessor) : base (name, start, end)
            {
                LeftNucsToTrim = leftNucsToTrim;
                //RightNucsToTrim = rightNucsToTrim;
                _sequenceAccessor = sequenceAccessor;
            }

            /*public string GetWholeSequence()
            {
                var halfLength = Length/2;
                var leftoverNucs = Length - (halfLength)*2;
                return string.Format("{0}{1}", GetSubSequence(halfLength + leftoverNucs, true), 
                    GetSubSequence(halfLength, false));                
            }*/

            protected override string GetSubSequence(uint start, uint end, uint length, 
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)//uint length, bool fromLeft)
            {
                var intLength = (int) length;
                //var sb = new StringBuilder(intLength + RightNucsToTrim);
                var returnChars = new char[intLength];
                var paddedStart = start + LeftNucsToTrim - 1;
                var bytePosition = paddedStart / NucsPerByte; // calculate the starting position to read from.
                var leftSubStringIndex = (int)paddedStart % NucsPerByte; // calculate the leftovers
                var lastByte = (end + LeftNucsToTrim - 1) / NucsPerByte;

                var charPosition = 0;
                var tetNuc = ByteToTetraNucleotide(_sequenceAccessor.ReadByte(bytePosition++));
                for (var i = leftSubStringIndex; i < NucsPerByte && charPosition < returnChars.Length; i++)
                    returnChars[charPosition++] = tetNuc[i]; // copy the leftmost side
                //sb.Append(ByteToTetraNucleotide(_sequenceAccessor.ReadByte(bytePosition++)), leftSubStringIndex, 4 - leftSubStringIndex));
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

                for (; bytePosition <= lastByte; bytePosition++)
                {
                    tetNuc = ByteToTetraNucleotide(_sequenceAccessor.ReadByte(bytePosition++));
                    for (var i = 0; i < NucsPerByte && charPosition < returnChars.Length; i++) // charPos < returnChars.Length short circuits
                        returnChars[charPosition++] = tetNuc[i];
                } 
                //sb.Append(ByteToTetraNucleotide(_sequenceAccessor.ReadByte(bytePosition)));}
                //sb.Remove(intLength, sb.Length-intLength); //trim end.
                /*if (fromLeft)
                {
                    if (nucsLeftover > 0)
                        sb.Append(ByteToTetraNucleotide(_sequenceAccessor.ReadByte(lastByte)).Substring(0, 4 - nucsLeftover));
                }
                else if (RightNucsToTrim > 0)
                    sb.Append(ByteToTetraNucleotide(_sequenceAccessor.ReadByte(lastByte)).Substring(0, 4 - RightNucsToTrim));*/
                //return sb.ToString();   
                return new string(returnChars);
            }

            public void Dispose()
            {
                _sequenceAccessor.Dispose();
            }
        }

        private class MaskedSubcontig : NormalSubcontig
        {
            public MaskedSubcontig(string name, uint start, uint end, ushort leftNucsToTrim, //ushort rightNucsToTrim, 
                MemoryMappedViewAccessor sequenceAccessor)
                : base(name, start, end, leftNucsToTrim, sequenceAccessor) { }

            protected override string GetSubSequence(uint start, uint end, uint length, 
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                if (skipMasks) return String.Empty;
                var tempstr = base.GetSubSequence(start, end, length, ignoreMasks, skipMasks, skipNs);
 	            return ignoreMasks ? tempstr : tempstr.ToLower();
            }
        }

        private class NSubcontig : TwoBitGenomeSubcontig
        {
            public NSubcontig(string name, uint start, uint end) : base(name, start, end) { }

            protected override string GetSubSequence(uint start, uint end, uint length,
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                return skipNs ? "" : new string('N', (int)length);
            }
        }

        private class TwoBitGenomeContig : IGenomeContig
        {
            public string Name { get; private set; }
            public uint Length { get; private set; }
            public bool ContainsNs { get; private set; }
            public bool ContainsMaskedSequences { get; private set; }

            private readonly List<ITwoBitGenomeSubcontig> _subSequences;

            public TwoBitGenomeContig(string name, uint length, List<ITwoBitGenomeSubcontig> subSequences)
            {
                Name = name;
                Length = length;
                _subSequences = subSequences;
                ContainsNs = _subSequences.Any(s => s.GetType() == typeof(NSubcontig));
                ContainsMaskedSequences = _subSequences.Any(s => s.GetType() == typeof(MaskedSubcontig));
            }

            public string GetSequence(uint start, uint end,
                bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                if (start < 1 || end > Length) throw new OutOfContigRangeException(Name, start, end,
                    string.Format("Inside {0}'s GetSequence method.", GetType().Name));
                if (end < start) throw new ArgumentException("Sequence end was less than start, what's up with that?");
                int position, endingIndex;
                if (Length - end > start - 1) // end is closer to the middle, so do the first search using end
                {
                    endingIndex = BinarySearchForIndex(end, 0, _subSequences.Count - 1);
                    if (start > _subSequences[endingIndex].Start) // all contained in one subSequence;
                        return _subSequences[endingIndex].GetSequence(start, end, ignoreMasks, skipMasks, skipNs);
                    position = BinarySearchForIndex(start, 0, endingIndex - 1);
                }
                else
                {
                    position = BinarySearchForIndex(start, 0, _subSequences.Count - 1);
                    if (end < _subSequences[position].End) // all contained in one subSequence;
                        return _subSequences[position].GetSequence(start, end, ignoreMasks, skipMasks, skipNs);
                    endingIndex = BinarySearchForIndex(end, position + 1, _subSequences.Count - 1);
                }
                var returnChars = new char[end - start + 1];
                var charPosition = 0;
                var seqString = _subSequences[position].GetSequence(start, _subSequences[position++].End, ignoreMasks, skipMasks, skipNs);
                for (var i = 0; i < seqString.Length && charPosition < returnChars.Length; i++)
                    returnChars[charPosition++] = seqString[i];
                for (; position <= endingIndex; position++)
                {
                    seqString = _subSequences[position].GetSequence(_subSequences[position].Start, 
                        _subSequences[position].End, ignoreMasks, skipMasks, skipNs);
                    for (var i = 0; i < seqString.Length && charPosition < returnChars.Length; i++)
                        returnChars[charPosition++] = seqString[i];
                }
                //seqString += _subSequences[position].GetSequence(_subSequences[position].Start, end, ignoreMasks, ignoreNs);
                //return seqString;
                return new string(returnChars);
            }

            private int BinarySearchForIndex(uint position, int begin, int end)
            {
                if (begin == end) return begin; // already found it.
                var middle = (begin + end) / 2;
                if (position < _subSequences[middle].Start) // before middle
                    return BinarySearchForIndex(position, begin, middle - 1);
                if (position > _subSequences[middle].End) // after middle
                    return BinarySearchForIndex(position, middle + 1, end);
                return middle;
            }

            public char GetNucleotideAt(uint position)
            {
                return GetSequence(position, position)[0];
            }

            public IEnumerable<char> GetNucleotides(uint start, uint end, bool ignoreMasks = false, bool skipMasks = false, bool skipNs = false)
            {
                return GetSequence(start, end, ignoreMasks, skipMasks, skipNs);
            }
        }

        public bool IsZeroBasedCoordinates { get { return false; } }
        public bool SupportsMasking { get { return true; } }
        public bool SupportsNs { get { return true; } }
        public bool SupportsIupacAmbiguityCodes { get { return false; } }

        private readonly List<IGenomeContig> _contigs;
        public List<IGenomeContig> Contigs { get { return _contigs; } }

        private readonly MemoryMappedFile _mmf;
        private readonly bool _reverse = false;
        private string _twoBitFile = null;
        private readonly uint _reserved;
        public uint Reserved { get { return _reserved; }}
        private string _fileError { get { return _twoBitFile ?? _mmf.ToString(); } }
        private const uint ForwardSignature = 0x1A412743;
        private const uint ReverseSignature = 0x4327411A;
        private const uint TwoBitVersion = 0;

        internal static readonly byte[] EndiannessIndexedByteArray = new byte[256];
        
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
            int numSequences;
            List<Tuple<string, uint>> contigInfos;
            using (var stream = File.OpenRead(file))
            {
                using (var br = new BinaryReader(stream))
                {
                    ReadSignature(br);
                    ReadVersion(br);
                    numSequences = (int) ReadUint(br);
                    unchecked
                    {
                        _contigs = new List<IGenomeContig>(numSequences);
                    }
                    _reserved = ReadUint(br);
                    contigInfos = new List<Tuple<string, uint>>(numSequences);
                    for (var i = 0; i < numSequences; i++)
                        contigInfos.Add(ReadSequenceIndex(br));
                }
            }
            _mmf = MemoryMappedFile.CreateFromFile(file, FileMode.Open);

        }

        public TwoBitGenomeReader(FileInfo file) : this(file.FullName) { }

        private bool ReadSignature(BinaryReader br)
        {
            var signature = ReadUint(br, true); // only time it will need to be read literally.
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
            for (var i = 0; i < n; i++)
                byteArray[i] = EndiannessIndexedByteArray[br.ReadByte()];
            return byteArray;
        }

        private uint ReadUint(BinaryReader br, bool readLiteral = false)
        {
            return readLiteral ? br.ReadUInt32() : BitConverter.ToUInt32(ReadNBytes(br, 4), 0);
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
