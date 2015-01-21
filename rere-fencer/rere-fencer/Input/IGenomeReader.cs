﻿using System;
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
        List<IGenomeContig> Contigs { get; }
    }

    public interface IGenomeContig
    {
        string Name { get; }
        ulong Length { get; }
        bool ContainsNs { get; }
        bool ContainsMaskedSequences { get; }
        string GetSequence(ulong start, ulong end, bool ignoreMasks = false, bool ignoreNs = false);
        char GetNucleotideAt(ulong position);
        IEnumerable<char> GetNucleotides(ulong start, ulong end, bool ignoreMasks = false, bool ignoreNs = false); 
        //uint? FirstPositionOfSequence(string sequence, uint offset = 0, bool strict = true);
        //uint? LastPositionOfSequence(string sequence, uint offset = 0, bool strict = true);
    }
}
