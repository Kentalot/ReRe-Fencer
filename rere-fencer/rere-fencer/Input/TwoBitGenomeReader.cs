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


    public class TwoBitGenomeContig : IGenomeContig
    {
        public string Name { get { return _name; } } private readonly string _name;
        public uint Length { get { return _length; } } private readonly uint _length;
        public bool ContainsNs { get { return _containsNs; } } private readonly bool _containsNs;
        public bool ContainsMaskedSequences { get { return _containsMaskedSequences; } } private readonly bool _containsMaskedSequences;

        public uint Reserved { get { return _reserved; } } private readonly uint _reserved;
        private readonly IntervalTree<bool, uint> _specialBlocks;
        private bool _isLastReadInNBlock = false; // used for slight optimization.

        //public uint NBlockCount { get { return _nBlockCount; } }
        private readonly MemoryMappedViewAccessor _contigSequence;

        public TwoBitGenomeContig(uint length, uint[] nBlockStarts, uint[] nBlockSizes, 
            uint[] maskBlockStarts, uint[] maskBlockSizes, uint reserved,
            MemoryMappedViewAccessor contigSequence)
        {
            _length = length;
            _specialBlocks = new IntervalTree<bool, uint>();
            UpdateBlocksTree(_specialBlocks, nBlockStarts, nBlockSizes, true);
            UpdateBlocksTree(_specialBlocks, maskBlockSizes, maskBlockStarts, false);
            _reserved = reserved;
            _contigSequence = contigSequence;
        }

        private void UpdateBlocksTree(IntervalTree<bool, uint> treeToUpdate, uint[] blockStarts, uint[] blockSizes, bool isNBlock)
        {
            var length = blockStarts.Length;
            if (length != blockSizes.Length) 
                throw new InvalidTwoBitContigFormatException(Name, 
                    string.Format("The number of {0}BlockStarts and {0}BlockSizes differ!", isNBlock ? "n" : "mask"));
            for (var i = 0; i < length; i++)
                treeToUpdate.AddInterval(blockStarts[i], blockStarts[i] + blockSizes[i] - 1, isNBlock);
        }

        public string GetAByteOfSequence(uint start)
        {
            var abyte = TwoBitGenomeReader.EndiannessIndexedByteArray[_contigSequence.ReadByte(start)];
            var sequence = ((TwoBitGenomeReader.TetraNucleotide) abyte).ToString();
            var overlappingIntervals = _specialBlocks.GetIntervals(start, start + 3);
            if (overlappingIntervals.Count < 1) return sequence;
            var seqArray = sequence.ToCharArray();
            for (uint i = 0; i < 4; i++)
            {
                Interval<bool, uint> overlapping = null;
                try
                {
                    overlapping = overlappingIntervals.First(v => v.ContainsWithStartEnd(start + i));
                } catch (InvalidOperationException){}
                if (overlapping == null) continue;
                seqArray[i] = overlapping.Data ? 'N' : Char.ToLower(seqArray[i]);
            }
        }

        public string GetSequenceAt(uint start, uint end, bool ignoreNs = false)
        {
            throw new NotImplementedException();
        }

        public char GetNucleotideAt(uint position)
        {
            if (OverlapsNBlock(position)) return 'N';
            var nucByte = GetAByteOfSequence(position);
            nucByte = OverlapsMaskBlock(position) ? nucByte.ToLower() : nucByte;
            return nucByte[0];
        }

        public uint? FirstPositionOfSequence(string sequence, uint offset = 0, bool strict = true)
        {
            throw new NotImplementedException();
        }

        public uint? LastPositionOfSequence(string sequence, uint offset = 0, bool strict = true)
        {
            throw new NotImplementedException();
        }

        private bool OverlapsMaskBlock(uint start, int size)
        {
            return 
        }

        private bool OverlapsNBlock(uint position, int size)
        {
            throw new NotImplementedException();
        }
    }

    public class TwoBitGenomeReader : IGenomeReader, IDisposable
    {
        public bool SupportsMasking { get { return true; } }
        public bool SupportsNs { get { return true; } }
        public bool SupportsIupacAmbiguityCodes { get { return false; } }

        private readonly SortedList<string, IGenomeContig> _contigs;
        public SortedList<string, IGenomeContig> Contigs { get; private set; }

        private readonly MemoryMappedFile _mmf;
        private readonly SequenceInfo[] _sequenceInfos;
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
            TTTT = 0x00,TTTC,TTTA,TTTG,TTCT,TTCC,TTCA,TTCG,TTAT,TTAC,TTAA,TTAG,TTGT,TTGC,TTGA,TTGG,
            TCTT,TCTC,TCTA,TCTG,TCCT,TCCC,TCCA,TCCG,TCAT,TCAC,TCAA,TCAG,TCGT,TCGC,TCGA,TCGG,
            TATT,TATC,TATA,TATG,TACT,TACC,TACA,TACG,TAAT,TAAC,TAAA,TAAG,TAGT,TAGC,TAGA,TAGG,
            TGTT,TGTC,TGTA,TGTG,TGCT,TGCC,TGCA,TGCG,TGAT,TGAC,TGAA,TGAG,TGGT,TGGC,TGGA,TGGG,
            CTTT,CTTC,CTTA,CTTG,CTCT,CTCC,CTCA,CTCG,CTAT,CTAC,CTAA,CTAG,CTGT,CTGC,CTGA,CTGG,
            CCTT,CCTC,CCTA,CCTG,CCCT,CCCC,CCCA,CCCG,CCAT,CCAC,CCAA,CCAG,CCGT,CCGC,CCGA,CCGG,
            CATT,CATC,CATA,CATG,CACT,CACC,CACA,CACG,CAAT,CAAC,CAAA,CAAG,CAGT,CAGC,CAGA,CAGG,
            CGTT,CGTC,CGTA,CGTG,CGCT,CGCC,CGCA,CGCG,CGAT,CGAC,CGAA,CGAG,CGGT,CGGC,CGGA,CGGG,
            ATTT,ATTC,ATTA,ATTG,ATCT,ATCC,ATCA,ATCG,ATAT,ATAC,ATAA,ATAG,ATGT,ATGC,ATGA,ATGG,
            ACTT,ACTC,ACTA,ACTG,ACCT,ACCC,ACCA,ACCG,ACAT,ACAC,ACAA,ACAG,ACGT,ACGC,ACGA,ACGG,
            AATT,AATC,AATA,AATG,AACT,AACC,AACA,AACG,AAAT,AAAC,AAAA,AAAG,AAGT,AAGC,AAGA,AAGG,
            AGTT,AGTC,AGTA,AGTG,AGCT,AGCC,AGCA,AGCG,AGAT,AGAC,AGAA,AGAG,AGGT,AGGC,AGGA,AGGG,
            GTTT,GTTC,GTTA,GTTG,GTCT,GTCC,GTCA,GTCG,GTAT,GTAC,GTAA,GTAG,GTGT,GTGC,GTGA,GTGG,
            GCTT,GCTC,GCTA,GCTG,GCCT,GCCC,GCCA,GCCG,GCAT,GCAC,GCAA,GCAG,GCGT,GCGC,GCGA,GCGG,
            GATT,GATC,GATA,GATG,GACT,GACC,GACA,GACG,GAAT,GAAC,GAAA,GAAG,GAGT,GAGC,GAGA,GAGG,
            GGTT,GGTC,GGTA,GGTG,GGCT,GGCC,GGCA,GGCG,GGAT,GGAC,GGAA,GGAG,GGGT,GGGC,GGGA,GGGG
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
