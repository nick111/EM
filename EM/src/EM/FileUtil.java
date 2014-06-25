package EM;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class FileUtil {
	
	public static List<String[]> getStringsByQuot(String fileAddress){

		List<String[]> strList = new ArrayList<String[]>();

		  File file = new File(fileAddress);
		  try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String oneLineStr;
			while((oneLineStr = br.readLine()) != null){
				strList.add(oneLineStr.split(","));
				}

		  } catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return strList;
		
	}

	
	public static List<double[]> getDoublesByQuot(String fileAddress){
		List<String[]> strList = getStringsByQuot(fileAddress);
		List<double[]> doubleList = new ArrayList<double[]>();
		for(int i = 0; i< strList.size(); i++){
			double[] dblList = new double[strList.get(i).length];
			for(int k = 0; k < strList.get(i).length; k++){
				dblList[k] =  Double.parseDouble(strList.get(i)[k]);
			}
			doubleList.add(dblList);
		}
		return doubleList;
	}
	
	public static List<String[]> getStringsByTab(String fileAddress){

		List<String[]> strList = new ArrayList<String[]>();

		  File file = new File(fileAddress);
		  try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String oneLineStr;
			while((oneLineStr = br.readLine()) != null){
				strList.add(oneLineStr.split("\t"));
				}

		  } catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return strList;
		
	}

	
	public static List<double[]> getDoublesByTab(String fileAddress){
		List<String[]> strList = getStringsByTab(fileAddress);
		List<double[]> doubleList = new ArrayList<double[]>();
		for(int i = 0; i< strList.size(); i++){
			double[] dblList = new double[strList.get(i).length];
			for(int k = 0; k < strList.get(i).length; k++){
				dblList[k] =  Double.parseDouble(strList.get(i)[k]);
			}
			doubleList.add(dblList);
		}
		return doubleList;
	}
	
	
	public static void makeFileByQuot(List<String[]> strlstList, String fileAddress){
		  File file = new File(fileAddress);
		  
		  try {
			  if (!checkBeforeWritefile(file)){
				  BufferedWriter bw = new BufferedWriter(new FileWriter(file));
				  
				  for(int i = 0; i < strlstList.size(); i++){
					  String str = "";
					  for(int k = 0; k < strlstList.get(i).length; k++){
						  str = str + strlstList.get(i)[k] +",";
					  }
					  str = str.substring(0, str.length() - 1);
					  bw.write(str);
					  bw.newLine();
				  }
				  bw.close();
				  
				  
			}
		  }catch (IOException e) {
			  e.printStackTrace();
			  }
		  }
	
	public static void makeFileBySpace(List<String[]> strlstList, String fileAddress){
		  File file = new File(fileAddress);
		  
		  try {
			  if (!checkBeforeWritefile(file)){
				  BufferedWriter bw = new BufferedWriter(new FileWriter(file));
				  
				  for(int i = 0; i < strlstList.size(); i++){
					  String str = "";
					  for(int k = 0; k < strlstList.get(i).length; k++){
						  str = str + strlstList.get(i)[k] +" ";
					  }
					  str = str.substring(0, str.length() - 1);
					  bw.write(str);
					  bw.newLine();
				  }
				  bw.close();
				  
				  
			}
		  }catch (IOException e) {
			  e.printStackTrace();
			  }
		  }

	
	private static boolean checkBeforeWritefile(File file){
		if (file.exists()){
			if (file.isFile() && file.canWrite()){
				return true;
				}
			}
		return false;
		}
	
}
