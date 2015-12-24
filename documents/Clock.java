
 
/** 

@author Per Reedtz Thomsen. Copyright (C) 1995 Per Reedtz Thomsen

Based on Rachel Gollub's Original clock Applet.

This program tries to demonstrate the use of parameters in the 
Beta interface for Java, while also attempting to be useful.

The applet displays a digital clock. The display can be configured:
Font Family, size, and weight, as well as background and foreground 
colors can be controlled from the html page that calls the applet.

The format of the displayed string can be controlled as well, using the 
formatting fields from the UNIX date(1) command, with the following 
exceptions:

  j, u, W, n, t conversions are not supported.


Example of how to call the applet:


<applet code="Clock.class" width=170 height=50>
<param name="FontFamily" value="Helvetica">
<param name="FontSize" value="14">
<param name="FontWeight" value="bold">
<param name="BGCol" value="f0f0f0">
<param name="FGCol" value="36648b">
<param name="DateFmt" value="%A, %B %I">
</applet>

LEGALESE:
Permission to use, copy, modify, and distribute this software
and its documentation for NON-COMMERCIAL or COMMERCIAL purposes and
without fee is hereby granted. 

No express promises are made, as to the usability of this software 
for any purpose whatsoever. If you depend on this for anything, and 
it breaks, and you lose money, don't come crying to me... 

*/


import java.util.*;
import java.awt.*;
import java.applet.*;
// import java.io.*;

public class Clock extends Applet implements Runnable {
  // Globals
  Thread timer = null; // The timer thread
  Color foregroundCol; // The foreground color
  Color DefaultforegroundCol = new Color(0, 0, 0); // The default Foreground (Black)
  Color backgroundCol; // The background color
  Color DefaultbackgroundCol = new Color(192, 192, 192); // 'Browser Gray'
  Date dummy = new Date(); // The first date to format the output from
  Font showfont; // The font used for display
  String lastdate = dummy.toLocaleString(); // The date to erase
  String DateFmt = "%a, %B %e %T";
  String[] weekDays = new String[7];
  String[] months = new String[12];  
  int earnedCredit = 1;  // The default earned credit.

/** 
Get a color parameter from the html page (<param ....>)
If the parameter is not given in the <applet> container,
return the default color.
*/
private int getColorFromParam(String parmName, Color defaultCol)
{
  Integer IntColor;
  String parmValue = this.getParameter(parmName); // Get the value
  if (parmValue != null) // Was the <param name="[parmName]" ...> there?
  { 
    Integer intColor = Integer.valueOf(parmValue, 16);
    return(intColor.intValue()); // return the 32-bit color value
  }
  return (defaultCol.getRGB()); // return the default color value
}


/** 
Prepare everything, so we will run correctly.
This includes getting and setting the parameters
obtained from the <param ...> tags.
*/ 
public void init()
{
  String FontFamily = "TimesRoman"; // The default Font Family
  int FontSize = 12; // The default Font Size
  int FontWeight = Font.PLAIN; // The default Font Weight
  // Set up the colors to use (FG & BG)
  backgroundCol = new Color(this.getColorFromParam("BGCol", DefaultbackgroundCol));
  foregroundCol = new Color(this.getColorFromParam("FGCol", DefaultforegroundCol));
  this.setBackground(backgroundCol); // Set the background Color

  // If the "FontFamily" param is present
  if (this.getParameter("FontFamily") != null)
  {
    FontFamily = this.getParameter("FontFamily"); // Change the Font Family
  }

  // If the "FontSize" param is present
  if (this.getParameter("FontSize") != null)
  {
    // Convert the value string into an Integer
    Integer IntFontSize = new Integer(this.getParameter("FontSize"));
    // and then to an int
    FontSize = IntFontSize.intValue();
  }

  // Tom's stuff
  // If the "earnedCredit" param is present
  if (this.getParameter("earnedCredit") != null)
  {
    // Convert the value string into an Integer
    Integer IntEarnedCredit = new Integer(this.getParameter("earnedCredit"));
    // and then to an int
    earnedCredit = IntEarnedCredit.intValue();
  }

  // If the "FontWeight" param is present
  if (this.getParameter("FontWeight") != null)
  {
    String SFontWeight = new String(this.getParameter("FontWeight"));
    SFontWeight = SFontWeight.toLowerCase(); // convert to lowercase
    // Set the appropriate values. NOTE: You can't or values together...
    if (SFontWeight.equals("bold")) { FontWeight = Font.BOLD; }
    if (SFontWeight.equals("italic")) { FontWeight = Font.ITALIC; }
    if (SFontWeight.equals("plain")) { FontWeight = Font.PLAIN; }
  } 
  showfont = new Font(FontFamily, FontWeight, FontSize); // Create the desired font.

  // if the "DateFmt" param is present
  if (this.getParameter("DateFmt") != null)
  {
    DateFmt = this.getParameter("DateFmt");
  }
  
  // Initialize the Formatting Arrays for date display
  
  // Days (Abbreviated is obtained by x.substring(0,3))
  weekDays[0] = new String("Sunday");
  weekDays[1] = new String("Monday");
  weekDays[2] = new String("Tuesday");
  weekDays[3] = new String("Wednesday");
  weekDays[4] = new String("Thursday");
  weekDays[5] = new String("Friday");
  weekDays[6] = new String("Saturday");

  // Months (Abbreviated is obtained by x.substring(0,3))
  months[0] = new String("January");
  months[1] = new String("February");
  months[2] = new String("March");
  months[3] = new String("April");
  months[4] = new String("May");
  months[5] = new String("June");
  months[6] = new String("July");
  months[7] = new String("August");
  months[8] = new String("September");
  months[9] = new String("October");
  months[10] = new String("November");
  months[11] = new String("December");

}


// Left Pad a number
private String padElement(int expr, char padChar)
{
  String result = "";
  // I'm just padding 2 digit numbers
  if (expr < 10) result = result.concat(String.valueOf(padChar));
  result = result.concat(String.valueOf(expr));
  return(result);
}

// Format a date according to the formatting string.
private String formatDate(String fmt, Date d)
{
  String formattedDate = "Projected Completion Date: "; // Start with an empty string
  
  // Retrieve the specific date information
  int year = d.getYear();
  int longYear = year + 1900;
  int month = d.getMonth();
  int monthDay = d.getDate();
  int weekDay = d.getDay();
  int hour = d.getHours();
  int minute = d.getMinutes();
  int second = d.getSeconds();

  // Adjust the year, if after 2000
  year = year > 99 ? year - 100 : year;
  int US_Hour = hour < 13 ? hour : hour - 12;

  // Loop through the format string
  for(int i = 0; i < fmt.length(); i++)
  {
    if (fmt.charAt(i) == '%') // We've hit a formatting command...
    {
      i++; // Move past the '%' sign
      if (fmt.length() <= i) // If the last char of the format is a lone '%'
      {
        formattedDate = formattedDate.concat("?");
        continue; // Forget the rest of the loop
      }
      // Figure out the format.
      switch (fmt.charAt(i))
      {
      case 'a': // Short Weekday name
        formattedDate = formattedDate.concat(weekDays[weekDay].substring(0,3));
        break;
      case 'A': // Long Weekday name
        formattedDate = formattedDate.concat(weekDays[weekDay]);
        break;
      case 'b': // Short Month name
      case 'h': // It's alias
        formattedDate = formattedDate.concat(months[month].substring(0,3));
        break;
      case 'B': // Long Month name
        formattedDate = formattedDate.concat(months[month]);
        break;
      case 'c': // The locale time/date string
        formattedDate = formattedDate.concat(d.toLocaleString());
        break;
      case 'C': // The default time/date string
        formattedDate = formattedDate.concat(d.toString());
        break;
      case 'd': // 2 digit month number
        formattedDate = formattedDate.concat(padElement(monthDay, '0'));
        break;
      case 'D': // Shortcut for %m/%d/%y
        formattedDate = formattedDate.concat(padElement(month + 1, '0'));
        formattedDate = formattedDate.concat(String.valueOf('/'));
        formattedDate = formattedDate.concat(padElement(monthDay, '0'));
        formattedDate = formattedDate.concat(String.valueOf('/'));
        formattedDate = formattedDate.concat(padElement(year, '0'));
        break;
      case 'e': // Month Number
        formattedDate = formattedDate.concat(padElement(monthDay, ' '));
        break;
      case 'H': // Hour -- 00 to 23
        formattedDate = formattedDate.concat(padElement(hour, '0'));
        break;
      case 'I': // Hour -- 01 to 12
        formattedDate = formattedDate.concat(padElement(US_Hour, '0'));
        break;
	// case 'j': // day of year 001 to 366 ; left out 
        // (java doesn't have the d.getYearDay() function.
      case 'm': // Month numbers -- 01 to 12
        formattedDate = formattedDate.concat(padElement(month + 1, '0'));
        break;
      case 'M': // Minutes -- 00 to 59
        formattedDate = formattedDate.concat(padElement(minute, '0'));
        break;
      // case 'n': // Insert a newline; I guess drawString doesn't do \x stuff
      //   formattedDate = formattedDate.concat("\n");
      //   break;
      case 'p': // AM or PM
        formattedDate = formattedDate.concat(String.valueOf((hour < 12 ? "AM" : "PM")));
        break;
      case 'r': // Shortcut for %I:%M:%S %p
        formattedDate = formattedDate.concat(padElement(US_Hour, '0'));
        formattedDate = formattedDate.concat(String.valueOf(':'));
        formattedDate = formattedDate.concat(padElement(minute, '0'));
        formattedDate = formattedDate.concat(String.valueOf(':'));
        formattedDate = formattedDate.concat(padElement(second, '0'));
        formattedDate = formattedDate.concat(String.valueOf(' '));
        formattedDate = formattedDate.concat(String.valueOf((hour < 12 ? "AM" : "PM")));
        break;
      case 'R': // Shortcut for %H:%M
        formattedDate = formattedDate.concat(padElement(hour, '0'));
        formattedDate = formattedDate.concat(String.valueOf(':'));
        formattedDate = formattedDate.concat(padElement(minute, '0'));
        break; 
      case 'S': // Second -- 00 to 61 (leap seconds)
        formattedDate = formattedDate.concat(padElement(second, '0'));
        break;
      // case 't': // Insert a tab character; It seems drawString doesn't do \x stuff
      //   formattedDate = formattedDate.concat("\t");
      //   break;
      case 'T': // Shortcut for %H:%M:%S
        formattedDate = formattedDate.concat(padElement(hour, '0'));
        formattedDate = formattedDate.concat(String.valueOf(':'));
        formattedDate = formattedDate.concat(padElement(minute, '0'));
        formattedDate = formattedDate.concat(String.valueOf(':'));
        formattedDate = formattedDate.concat(padElement(second, '0'));
        break; 
	// 'u' and 'W' not supported.
      case 'w': // day of week in numbers; Sunday = 0   
        formattedDate = formattedDate.concat(String.valueOf(weekDay));
        break;
	// 'x' and 'X' not supported
      case 'y': // short year (00 - 99)
        formattedDate = formattedDate.concat(padElement(year, '0'));
        break;
      case 'Y': // long year (1996)
        formattedDate = formattedDate.concat(padElement(longYear, '0'));
        break;
	// 'Z' not supported
      case '%': // Just in case you want to show the '%' sign
        formattedDate = formattedDate.concat("%");
        break;
      default:
        formattedDate = formattedDate.concat("??");
        break;
      }
    }
    else // A regular character
    {
      formattedDate = formattedDate.concat(String.valueOf(fmt.charAt(i)));
    }
  } // end for

  return(formattedDate);
}

/**
Paint is the main part of the program
First get the current date, and set the correct font;
then erase the old time/date, and draw the new one
*/
public void paint(Graphics g)
{
  String today;
  Date dat = new Date();
  // Tom Hales's modifications of display date.
  // long starting = 874120205000L;
  long starting =    876667709000L;

  long current = dat.getTime();
  double totalCredit = 27000. + 1.;  // 
  dat = new Date(starting + 
	((long)(totalCredit/((double) earnedCredit)))*(current-starting));

  g.setFont(showfont);

  today = formatDate(DateFmt, dat);
  // today = dat.toLocaleString();
  
  // Clean up the old text
  g.setColor(backgroundCol);
  g.drawString(lastdate, 5, 20);
  // Draw the new text
  g.setColor(foregroundCol);
  g.drawString(today, 5, 20); 
  // Get ready to erase the old text. 
  lastdate = today;
}

// start the thread
public void start()
{
  if(timer == null)
    {
      timer = new Thread(this);
      timer.start();
    }
}

/** 
When the timer Thread is set to null, the applet will stop;
see start()
*/
public void stop()
{
  timer = null;
}


public void run()
{

  // Sleep in the timer thread...
	
  //while (timer != null) {
    //try {timer.sleep(100);} catch (InterruptedException e){}
     repaint(); // and do the redraw.
  //}
   timer = null;
}


// Very simple update()
public void update(Graphics g)
{
  paint(g);
}
}
