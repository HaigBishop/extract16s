def analyze_alignment_file(file_path):
    with open(file_path, 'r') as f:
        # Skip the header line
        header = f.readline()
        
        # Read the rest of the file as the sequence
        sequence = ''.join(f.readlines())
        
        # Remove newlines
        sequence = sequence.replace('\n', '')
        
        # Count characters
        total_chars = len(sequence)
        
        # Count specific character types
        lowercase_count = sum(1 for c in sequence if c.islower())
        uppercase_count = sum(1 for c in sequence if c.isupper())
        dash_count = sequence.count('-')
        
        # Print results
        print(f"1. Num characters: {total_chars}")
        print(f"2. Num lowercases: {lowercase_count}")
        print(f"3. Num uppercases: {uppercase_count}")
        print(f"4. Num dashes: {dash_count}")
        
        # Extract V4 region (positions 531-783)
        # We need to count characters (both upper and lowercase, excluding dashes)
        nucleotide_count = 0
        v4_start_pos = None
        v4_end_pos = None
        
        for i, char in enumerate(sequence):
            if char != '-':  # Only count nucleotides, not gaps
                nucleotide_count += 1
                
                if nucleotide_count == 531 and v4_start_pos is None:
                    v4_start_pos = i
                    
                if nucleotide_count == 783 and v4_end_pos is None:
                    v4_end_pos = i
                    break
        
        # Extract the V4 subsequence
        if v4_start_pos is not None and v4_end_pos is not None:
            v4_region = sequence[v4_start_pos:v4_end_pos+1]
            
            # Format: uppercase and no dashes
            v4_region_formatted = v4_region.upper().replace('-', '')
            
            print(f"\nV4 region (positions 531-783):")
            print(v4_region_formatted)
        else:
            print("\nCould not extract V4 region - sequence may be too short")
        
        # Extract V3-V4 region (positions 355-781)
        # We need to count characters (both upper and lowercase, excluding dashes)
        nucleotide_count = 0
        v3v4_start_pos = None
        v3v4_end_pos = None
        
        for i, char in enumerate(sequence):
            if char != '-':  # Only count nucleotides, not gaps
                nucleotide_count += 1
                
                if nucleotide_count == 355 and v3v4_start_pos is None:
                    v3v4_start_pos = i
                    
                if nucleotide_count == 781 and v3v4_end_pos is None:
                    v3v4_end_pos = i
                    break
        
        # Extract the V3-V4 subsequence
        if v3v4_start_pos is not None and v3v4_end_pos is not None:
            v3v4_region = sequence[v3v4_start_pos:v3v4_end_pos+1]
            
            # Format: uppercase and no dashes
            v3v4_region_formatted = v3v4_region.upper().replace('-', '')
            
            print(f"\nV3-V4 region (positions 355-781):")
            print(v3v4_region_formatted)
        else:
            print("\nCould not extract V3-V4 region - sequence may be too short")
            
        # Extract V4 region based on alignment indices (524-776)
        # Only capital letters and dashes count as alignment positions
        alignment_count = 0
        v4_align_start_pos = None
        v4_align_end_pos = None
        
        for i, char in enumerate(sequence):
            if char.isupper() or char == '-':  # Count only uppercase letters and dashes
                alignment_count += 1
                
                if alignment_count == 524 and v4_align_start_pos is None:
                    v4_align_start_pos = i
                    
                if alignment_count == 776 and v4_align_end_pos is None:
                    v4_align_end_pos = i
                    break
        
        # Extract the V4 subsequence by alignment indices
        if v4_align_start_pos is not None and v4_align_end_pos is not None:
            v4_align_region = sequence[v4_align_start_pos:v4_align_end_pos+1]
            
            # Format: uppercase and no dashes
            v4_align_region_formatted = v4_align_region.upper().replace('-', '')
            
            print(f"\nV4 region (alignment indices 524-776):")
            print(v4_align_region_formatted)
        else:
            print("\nCould not extract V4 region by alignment indices - sequence may be too short")


if __name__ == "__main__":
    analyze_alignment_file("bac_a2m.fna")
