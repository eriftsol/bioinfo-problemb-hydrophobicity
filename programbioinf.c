#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WINDOW_SIZE 9		// Size of the sliding window
#define CUTOFF 0.6		// Hydrophobicity cutoff for identifying transmembrane segments

double Aminoacids (char);	//Prototype

int
main (int argc, char *argv[])
{
  if (argc != 2)
    {				//Case of wrong number of arguments
      printf ("Usage: %s protein_sequence\n", argv[0]);
      return 1;
    }

  char *protein_sequence = argv[1];	//pointer to protein array
  int sequence_length = strlen (protein_sequence);
  double avg_hydrophobicity = 0.0;
  int found_flag = 0;
  int digits = 0;


  // Declare an array of characters containing all  characters commonly found in protein sequences.
  char list[] = "aAcCdDeEfFgGhHiIkKlLmMnNpPqQrRsStTvVwWyY";

  // Initialize a variable to keep track of whether any characters were found
  // to not be in the list.
  int found = 1;

  // Iterate over the characters in the protein sequence.
  for (int i = 0; i < sequence_length; i++)
    {
      // Check if the current character is not in the list.
      if (!(strchr (list, protein_sequence[i])))
	{
	  // Print an error message if the character was not found.
	  printf
	    ("Character '%c' not found in list and it is in the position '%d'\n",
	     protein_sequence[i], i);
	  // Set the found variable to 0 to indicate that an invalid character was found.
	  found = 0;
	}
      // If the character was found in the list, continue to the next iteration.

      else
	{
	  continue;
	}
    }
  // If any invalid characters were found, print a message and exit the program.
  if (found == 0)
    {
      printf ("Stopping program because at least one character was invalid.\n");
      return 0;
    }



  for (int i = 0; i < sequence_length; i++)
    {				// Loop over all windows in the protein sequence
      digits = 0;
      avg_hydrophobicity = 0.0;

      for (int j = i; j < i + WINDOW_SIZE && j < sequence_length; j++)
	{			// Calculate the average hydrophobicity of the residues in the window
	  avg_hydrophobicity =
	    avg_hydrophobicity + Aminoacids (protein_sequence[j]);
	  digits++;
	}
      printf ("\n");
      avg_hydrophobicity /= WINDOW_SIZE;

      if ((avg_hydrophobicity > CUTOFF) && (digits == WINDOW_SIZE))
	{			//checking # of digits to avoid sample with less aminoacids than sliding window 
	  printf
	    ("Potential Transmembrane Segment found starting at position %d, with average hydrophobicity %f which is over the CUTOFF (%f) and its aminoacids are: ",
	     i, avg_hydrophobicity, CUTOFF);
	  for (int k = i; k < i + WINDOW_SIZE && k <= sequence_length; k++)
	    {			//prints of aminoacids
	      printf ("%c", protein_sequence[k]);
	    }
	  printf ("\n");
	  found_flag = 1;
	}
    }
  if (found_flag == 0)
    {
      printf ("No Potential Transmembrane Segment found!\n");
    }
  return 0;
}

double
Aminoacids (char ch)
{
  double x;
  switch (ch)
    {
    case 'A':
      x = 0.616;
      break;
    case 'a':
      x = 0.616;
      break;

    case 'C':
      x = 0.680;
      break;
    case 'c':
      x = 0.680;
      break;

    case 'D':
      x = 0.028;
      break;
    case 'd':
      x = 0.028;
      break;

    case 'E':
      x = 0.043;
      break;
    case 'e':
      x = 0.043;
      break;

    case 'F':
      x = 1;
      break;
    case 'f':
      x = 1;
      break;

    case 'G':
      x = 0.501;
      break;
    case 'g':
      x = 0.501;
      break;

    case 'H':
      x = 0.165;
      break;
    case 'h':
      x = 0.165;
      break;
    case 'I':
      x = 0.943;
      break;
    case 'i':
      x = 0.943;
      break;

    case 'K':
      x = 0.283;
      break;
    case 'k':
      x = 0.283;
      break;

    case 'L':
      x = 0.943;
      break;
    case 'l':
      x = 0.943;
      break;

    case 'M':
      x = 0.738;
      break;
    case 'm':
      x = 0.738;
      break;
    case 'N':
      x = 0.236;
      break;
    case 'n':
      x = 0.236;
      break;

    case 'P':
      x = 0.711;
      break;
    case 'p':
      x = 0.711;
      break;

    case 'Q':
      x = 0.251;
      break;
    case 'q':
      x = 0.251;
      break;

    case 'R':
      x = 0.000;
      break;
    case 'r':
      x = 0.000;
      break;
    case 'S':
      x = 0.359;
      break;
    case 's':
      x = 0.359;
      break;

    case 'T':
      x = 0.450;
      break;
    case 't':
      x = 0.450;
      break;

    case 'V':
      x = 0.825;
      break;
    case 'v':
      x = 0.825;
      break;

    case 'W':
      x = 0.878;
      break;
    case 'w':
      x = 0.878;
      break;

    case 'Y':
      x = 0.880;
      break;
    case 'y':
      x = 0.880;
      break;

    }
  return x;
}