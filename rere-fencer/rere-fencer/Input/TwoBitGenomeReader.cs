using System;
using System.CodeDom;
using System.Collections.Generic;
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
        private readonly MemoryMappedFile _mmf;
        private readonly SequenceInfo[] _sequenceInfos;
        private readonly bool _reverse = false;
        private string _twoBitFile = null;
        private readonly uint _reserved;
        public uint Reserved { get { return _reserved; }}
        private string _fileError { get { return _twoBitFile ?? _mmf.ToString(); } }
        internal static readonly byte[] MappedByteArray = new byte[256];

        private class SequenceInfo
        {
            public class SequenceData
            {
                private readonly uint _size;
                //public uint Size { get { return _size; } }
                private readonly uint _nBlockCount;
                private readonly uint[] _nBlockStarts;
                private readonly uint[] _nBlockSizes;
                private readonly uint _maskBlockCount;
                private readonly uint[] _maskBlockStarts;
                private readonly uint[] _maskBlockSizes;
                private readonly uint _reserved;
                public uint Reserved { get { return _reserved;} }

                //public uint NBlockCount { get { return _nBlockCount; } }
                private readonly MemoryMappedViewAccessor _contigSequence;

                public SequenceData(uint size, uint nBlockCount, 
                    MemoryMappedViewAccessor contigSequence)
                {
                    _size = size;
                    _nBlockCount = nBlockCount;
                    _contigSequence = contigSequence;
                }

                public string GetTetraNucleotideAt(long start)
                {
                    var abyte = MappedByteArray[_contigSequence.ReadByte(start)];
                    return ((TetraNucleotide) abyte).ToString();
                }

                public char GetBaseAt(long position)
                {
                    return GetTetraNucleotideAt(position)[0];
                }
            }
            private readonly string _sequenceName;
            public string SequenceName { get { return _sequenceName; } }
            private readonly uint _origin;
            public uint Origin { get { return _origin; } }
            public SequenceData Data { get; private set; }

            

            
        }

        //private static char[] bit_chars = {'T', 'C', 'A', 'G'};

        enum TetraNucleotide : byte
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
            using (var stream = File.OpenRead(file))
            {
                using (var br = new BinaryReader(stream))
                {
                    ReadSignature(br);
                    ReadVersion(br);
                    _sequenceInfos = new SequenceInfo[ReadUint(br)];
                    _reserved = ReadUint(br);
                }
            }

            _mmf = MemoryMappedFile.CreateFromFile(file, FileMode.Open);
        }

        public TwoBitGenomeReader(FileInfo file) : this(file.FullName) { }

        private bool ReadSignature(BinaryReader br)
        {
            var signature = ReadUint(br);
            bool reverse = false;
            if (signature == 0x1A412743) reverse = true;
            else if (signature != 0x4327411A) 
                throw new InvalidTwoBitFileException(_fileError, "Invalid Signature found: " + signature);
            for (byte i = 0; i < MappedByteArray.Length; i++)
                MappedByteArray[i] = reverse ? ReverseByte(i) : i;
            return reverse;
        }

        private void ReadVersion(BinaryReader br)
        {
            var version = ReadUint(br);
            if (version != 0) throw new InvalidTwoBitFileException(_fileError, "Version is something other than 0: " + version);
        }

        private byte[] ReadNBytes(BinaryReader br, int n)
        {
            var byteArray = new byte[n];
            for (var i = 0; i < n; i++)
                byteArray[i] = MappedByteArray[br.ReadByte()];
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
