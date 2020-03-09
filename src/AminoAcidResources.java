class AminoAcidResources{
  public static char getAminoAcidFromCodon(String s){
    //if this is not a 3 character string, return NULL
    if(s.length()!=3) return (char)0;
    s = s.toUpperCase();
    if(s.charAt(0) == 'A'){
      if(s.charAt(1) == 'A'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G') return 'K';
        if(s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'N';
      }
      if(s.charAt(1) == 'C'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G' || s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'T';
      }
      if(s.charAt(1) == 'G'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G') return 'R';
        if(s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'S';
      }
      if(s.charAt(1) == 'U'){
        if(s.charAt(2) == 'G') return 'M';
        if(s.charAt(2) == 'A' || s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'I';
      }
    }
    if(s.charAt(0) == 'C'){
      if(s.charAt(1) == 'A'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G') return 'Q';
        if(s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'H';
      }
      if(s.charAt(1) == 'C'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G' || s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'P';
      }
      if(s.charAt(1) == 'G'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G' || s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'R';
      }
      if(s.charAt(1) == 'U'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G' || s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'L';
      }
    }
    if(s.charAt(0) == 'G'){
      if(s.charAt(1) == 'A'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G') return 'E';
        if(s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'D';
      }
      if(s.charAt(1) == 'C'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G' || s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'A';
      }
      if(s.charAt(1) == 'G'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G' || s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'G';
      }
      if(s.charAt(1) == 'U'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G' || s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'V';
      }
    }
    if(s.charAt(0) == 'U'){
      if(s.charAt(1) == 'A'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G') return '*';//STOP
        if(s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'T';
      }
      if(s.charAt(1) == 'C'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G' || s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'S';
      }
      if(s.charAt(1) == 'G'){
        if(s.charAt(2) == 'A' ) return '*'; //STOP
        if(s.charAt(2) == 'G') return 'W';
        if(s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'C';
      }
      if(s.charAt(1) == 'U'){
        if(s.charAt(2) == 'A' || s.charAt(2) == 'G') return 'L';
        if(s.charAt(2) == 'C' || s.charAt(2) == 'U') return 'F';
      }
    }
    // If we ended up not finding all 3 characters, return NULL
    return (char)0;
  }
  public static String[] getCodonListForAminoAcid(char a){
    a = Character.toUpperCase(a);
    if(a == 'A') return new String[]{"GCG","GCA","GCC","GCU"};
    else if(a == 'C') return new String[]{"UGC","UGU"};
    else if(a == 'D') return new String[]{"GAC","GAU"};
    else if(a == 'E') return new String[]{"GAG","GAA"};
    else if(a == 'F') return new String[]{"UUC","UUU"};
    else if(a == 'G') return new String[]{"GGG","GGA","GGC","GGU"};
    else if(a == 'H') return new String[]{"CAC","CAU"};
    else if(a == 'I') return new String[]{"AUA","AUC","AUU"};
    else if(a == 'K') return new String[]{"AAG","AAA"};
    else if(a == 'L') return new String[]{"CUG","CUA","CUC","CUU","UUG","UUA"};
    else if(a == 'M') return new String[]{"AUG"};
    else if(a == 'N') return new String[]{"AAC","AAU"};
    else if(a == 'P') return new String[]{"CCG","CCA","CCC","CCU"};
    else if(a == 'Q') return new String[]{"CAG","CAA"};
    else if(a == 'R') return new String[]{"AGG","AGA","CGG","CGA","CGC","CGU"};
    else if(a == 'S') return new String[]{"AGC","AGU","UCG","UCA","UCC","UCU"};
    else if(a == 'T') return new String[]{"ACG","ACA","ACC","ACU","UAC","UAU"};
    else if(a == 'V') return new String[]{"GUG","GUA","GUC","GUU"};
    else if(a == 'W') return new String[]{"UGG"};
    // if character passed is not an amino acid, return an empty array
    return new String[]{};
  }
}