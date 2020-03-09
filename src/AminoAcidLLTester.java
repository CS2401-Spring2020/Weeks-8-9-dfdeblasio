import static org.junit.jupiter.api.Assertions.*;
import java.util.concurrent.TimeUnit;

import org.junit.jupiter.api.Test;
import java.util.Random;

public class AminoAcidLLTester {

  @Test
  public void smallTest() {
    AminoAcidLL list = AminoAcidLL.createFromRNASequence("AAAAAA");
    assertArrayEquals(new char[]{'K'},list.aminoAcidList());
    assertArrayEquals(new int[]{2},list.aminoAcidCounts());
  }

  @Test
  public void smallTest2() {
    AminoAcidLL list = AminoAcidLL.createFromRNASequence("CCCAAA");
    assertArrayEquals(new char[]{'P','K'},list.aminoAcidList());
    assertArrayEquals(new int[]{1,1},list.aminoAcidCounts());
  }

  @Test
  public void allPossibleAAs() {
    char[] aa = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W'};
    int[] counts = new int[aa.length];
    String s = "";
    for(int i=0; i<aa.length; i++){
      String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aa[i]);
      counts[i] = codonList.length;
      for(int l=0; l<codonList.length; l++){
        s += codonList[l];
      }
    }
    AminoAcidLL list = AminoAcidLL.createFromRNASequence(s);
    assertArrayEquals(aa,list.aminoAcidList());
    assertArrayEquals(counts,list.aminoAcidCounts());
  }

  @Test
  public void allPossibleAAsRev() {
    char[] aa = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W'};
    char[] aarev = new char[aa.length];
    int[] counts = new int[aa.length];
    String s = "";
    for(int i=aa.length-1; i>=0; i--){
      aarev[i] = aa[aa.length-1-i];
      String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aa[i]);
      counts[i] = codonList.length;
      for(int l=0; l<codonList.length; l++){
        s += codonList[l];
      }
    }
    AminoAcidLL list = AminoAcidLL.createFromRNASequence(s);
    assertFalse(list.isSorted());
    assertArrayEquals(aarev,list.aminoAcidList());
    list = AminoAcidLL.sort(list);
    assertTrue(list.isSorted());
    assertArrayEquals(aa,list.aminoAcidList());
    assertArrayEquals(counts,list.aminoAcidCounts());
  }
  
  @Test
  public void allPossibleAAsRevCompare() {
    char[] aa = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W'};
    String sR = "";
    for(int i=aa.length-1; i>=0; i--){
      String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aa[i]);
      for(int l=0; l<codonList.length; l++){
        sR += codonList[l];
      }
    }
    String s = "";
    for(int i=0; i<aa.length; i++){
      String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aa[i]);
      for(int l=0; l<codonList.length; l++){
        s += codonList[l];
      }
    }
    AminoAcidLL list = AminoAcidLL.createFromRNASequence(s);
    AminoAcidLL listRev = AminoAcidLL.createFromRNASequence(sR);
    assertTrue(list.isSorted());
    assertFalse(listRev.isSorted());
    listRev = AminoAcidLL.sort(listRev);
    assertTrue(listRev.isSorted());
    assertEquals(0,list.aminoAcidCompare(listRev));
    assertEquals(0,listRev.aminoAcidCompare(list));
    assertEquals(0,list.codonCompare(listRev));
    assertEquals(0,listRev.codonCompare(list));
  }
  
  @Test
  public void allPossibleAAsingle() {
    char[] aa = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W'};
    String sR = "";
    for(int i=aa.length-1; i>=0; i--){
      String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aa[i]);
      for(int l=0; l<codonList.length; l++){
        sR += codonList[l];
      }
    }
    String s = "AAA";
    AminoAcidLL list = AminoAcidLL.createFromRNASequence(s);
    AminoAcidLL listRev = AminoAcidLL.createFromRNASequence(sR);
    assertEquals(1,list.aminoAcidList().length);
    assertEquals(aa.length,listRev.aminoAcidList().length);
    assertTrue(list.isSorted());
    assertFalse(listRev.isSorted());
    listRev = AminoAcidLL.sort(listRev);
    assertTrue(listRev.isSorted());

    assertEquals(1,list.aminoAcidList().length);
    assertEquals(aa.length,listRev.aminoAcidList().length);
    
    assertEquals(60,list.aminoAcidCompare(listRev));
    assertEquals(60,listRev.aminoAcidCompare(list));
    assertEquals(60,list.codonCompare(listRev));
    assertEquals(60,listRev.codonCompare(list));
  }
  
  @Test
  public void allPossibleAARandomChanges() {
    String s1 = "UGGGUCACGUCGCGGCAACCUAAUAUGUUAAAAAUACACGGAUUUGAAGAUUGUGCA";
    String s2 = "UGGGUCUACUCACGUCAACCAAAUAUGCUCAAAAUUCAUGGCUUUGAGGACUGCGCC";
    AminoAcidLL list1 = AminoAcidLL.createFromRNASequence(s1);
    AminoAcidLL list2 = AminoAcidLL.createFromRNASequence(s2);
    assertEquals(19,list1.aminoAcidList().length);
    assertEquals(19,list2.aminoAcidList().length);
    assertFalse(list1.isSorted());
    assertFalse(list2.isSorted());
    list1 = AminoAcidLL.sort(list1);
    list2 = AminoAcidLL.sort(list2);
    assertTrue(list1.isSorted());
    assertTrue(list2.isSorted());
    
    assertEquals(0,list1.aminoAcidCompare(list2));
    assertEquals(0,list2.aminoAcidCompare(list1));
    assertNotEquals(0,list1.codonCompare(list2));
    assertNotEquals(0,list2.codonCompare(list1));
  }
}
