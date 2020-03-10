class AminoAcidLL{
  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;

  AminoAcidLL(){

  }


  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon 
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon){
    aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
    codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
    counts = new int[codons.length];
    next = null;
    incrementCodon(inCodon);
  }


  /********************************************************************************************/
  /* Given a codon string, finds the index for that codon and increments the count */
  private void incrementCodon(String inCodon){
    for(int i=0; i<codons.length; i++){
      if(codons[i].equals(inCodon.toUpperCase())){
        counts[i]++;
      }
    }
  }


  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops, 
   * if not passes the task to the next node. 
   * If there is no next node, add a new node to the list that would contain the codon. 
   */
  private void addCodon(String inCodon){
    if(AminoAcidResources.getAminoAcidFromCodon(inCodon) == aminoAcid){
      incrementCodon(inCodon);
    }else if(next != null){
      next.addCodon(inCodon);
    }else {
      next = new AminoAcidLL(inCodon);
    }
  }


  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
    int total = 0;
    for(int i=0; i<codons.length; i++){
      total += counts[i];
    }
    return total;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts. 
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){
    if(!isSorted() && inList.isSorted()) return -1;

    int diff = 0;
    //System.out.println(aminoAcid + " " + inList.aminoAcid);

    // if the comparison list is empty, the difference is the sum of the rest of the lists' counts
    if(inList == null){
      diff += totalCount();
      if(next != null){
        diff += next.aminoAcidCompare(inList);
      }


    }
    // the two nodes match, add the difference in counts from the rest of the list to this one and return the total. 
    else if(inList.aminoAcid == aminoAcid){
      diff = Math.abs(totalCount()-inList.totalCount());
      if(next != null){
        diff += next.aminoAcidCompare(inList.next);
      }if(next == null && inList.next != null){
        diff += aminoAcidCompare(inList.next);
      }

      // need to find out if something later in my list matches the current, 
      // but since I don't have a match my total count gets added to the difference
    }

    else if(next != null && aminoAcid < inList.aminoAcid){
      diff += totalCount();
      if(next != null){
        diff += next.aminoAcidCompare(inList);
      }

      // need to find out if something later in *their* list matches me, 
      // but since they don't have a match in my list their total count gets added to the difference
      // also if I don't have anything else keep adding their total to the difference
    }

    else if(next == null || aminoAcid > inList.aminoAcid){
      diff += inList.totalCount();
      if(inList.next != null){
        diff += aminoAcidCompare(inList.next);
      }

    }
    
    return diff;
  }

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }


  /********************************************************************************************/
  /* Same ad above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList){
    if(!isSorted() && inList.isSorted()) return -1;
    int diff = 0;

    // if the comparison list is empty, the difference is the sum of the rest of the lists' counts
    if(inList == null){
      diff += totalCount();
      if(next != null){
        diff += next.codonCompare(inList);
      }
    }

    // the two nodes match, add the difference in counts from the rest of the list to this one and return the total. 
    else if(aminoAcid == inList.aminoAcid){
      diff = codonDiff(inList);
      if(next != null){
        diff += next.codonCompare(inList.next);
      }if(next == null && inList.next != null){
        diff += codonCompare(inList.next);
      }


      // need to find out if something later in my list matches the current, 
      // but since I don't have a match my total count gets added to the difference
    }else if(next != null && aminoAcid < inList.aminoAcid){
      diff += totalCount();
      if(next != null){
        diff +=  next.codonCompare(inList);
      }

      // need to find out if something later in *their* list matches me, 
      // but since they don't have a match in my list their total count gets added to the difference
      // also if I don't have anything else keep adding their total to the difference
    }else if(next == null || aminoAcid > inList.aminoAcid){
      diff += inList.totalCount();
      if(inList.next != null){
        diff +=  codonCompare(inList.next);
      }
    }

    return diff;
  }


  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){
    if(next == null) {
      return new char[] {aminoAcid};
    }
    char[] temp = next.aminoAcidList();
    char[] rtn = new char[temp.length+1];
    for(int i=0;i<temp.length;i++) {
      rtn[i+1] = temp[i];
    }
    rtn[0] = aminoAcid;
    return rtn;

  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){
    if(next == null) {
      return new int[] {totalCount()};
    }
    int[] temp = next.aminoAcidCounts();
    int[] rtn = new int[temp.length+1];
    for(int i=0;i<temp.length;i++) {
      rtn[i+1] = temp[i];
    }
    rtn[0] = totalCount();
    return rtn;
  }


  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){
    boolean sorted = true;
    if(next != null){
      if(aminoAcid < next.aminoAcid){
        sorted &= next.isSorted();
      }else{
        sorted = false;
      }
    }
    return sorted;
  }


  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence){
    AminoAcidLL rtn = null;
    if(inSequence.length() >= 3 && AminoAcidResources.getAminoAcidFromCodon(inSequence.substring(0, 2)) != '*'){
      rtn = new AminoAcidLL(inSequence.substring(0, 3));
      boolean keepGoing = true;
      for(int i=3; i<inSequence.length()-2 && keepGoing; i+=3){
        if(AminoAcidResources.getAminoAcidFromCodon(inSequence.substring(i, i+3)) != '*'){
          rtn.addCodon(inSequence.substring(i, i+3));
        }else{
          keepGoing = false;
        }
      }
    }
    return rtn;
  }

  /********************************************************************************************/
  /* sorts a list by amino acid character 
   * using  either merge or insertion sort*/
  public static AminoAcidLL sort(AminoAcidLL inList){
	  return mergeSort(inList);
	  //return insertSort(inList);
  }

  /********************************************************************************************/
  /* sorts a list by amino acid character using merge sort*/
  public static AminoAcidLL mergeSort(AminoAcidLL inList){
    if(inList.next == null) return inList;
    
    // count the number of nodes 
    int ctr = 0;
    AminoAcidLL current = inList;
    while(current != null) {
    	current = current.next;
    	ctr++;
    }
    
    // create right to point to the second half (and clip the two lists).
    AminoAcidLL left = inList;
    current = inList;
    for(int i=0; i<ctr/2-1; i++) {
    	current = current.next;
    }
    AminoAcidLL right = current.next;
    current.next = null;
    
    // sort the shorter lists
    left = mergeSort(left);
    right = mergeSort(right);
    
    //merge 
    
    //set fist node
    AminoAcidLL rtn;
    if(left.aminoAcid < right.aminoAcid) {
    	rtn = left;
    	left = left.next;
    	rtn.next = null;
    }else {
    	rtn = right;
    	right = right.next;
    	rtn.next = null;
    }
    AminoAcidLL back = rtn;
    
    // combine the rest of the nodes
    while(left != null && right != null) {
    	if(left.aminoAcid < right.aminoAcid) {
    		back.next = left;
    		left = left.next;
    		back = back.next;
    		back.next = null;
    	}else {
    		back.next = right;
    		right = right.next;
    		back = back.next;
    		back.next = null;
    	}
    }
    //these cant both be true;
    if(left!=null) back.next = left;
    if(right!=null) back.next = right;
    
    return rtn;
  }
  
  /********************************************************************************************/
  /* sorts a list by amino acid character using indertion sort*/
  public static AminoAcidLL insertSort(AminoAcidLL inList){
    //using insertionSort
    if(inList == null) return inList;
    AminoAcidLL rtn = inList;
    AminoAcidLL curr = inList.next;

    rtn.next = null;
    while(curr != null) {
      //if cur is less than the current front add it as rtn
      if(rtn.aminoAcid > curr.aminoAcid) {
        AminoAcidLL temp = curr.next;
        curr.next = rtn;
        rtn = curr;
        curr = temp;
      }else{
        AminoAcidLL temp = rtn.next;
        AminoAcidLL prevtemp = rtn;
        while(temp != null && temp.aminoAcid < curr.aminoAcid) {
          prevtemp = temp;
          temp = temp.next;
        }
        if(temp == null) {
          prevtemp.next = curr;
          temp = curr.next;
          curr.next = null;
          curr = temp;
        }else {
          AminoAcidLL temp2 = curr.next;
          curr.next = temp;
          prevtemp.next = curr;
          curr = temp2;
        }
      }
    }
    return rtn;
  }
}