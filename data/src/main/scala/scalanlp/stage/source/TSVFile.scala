/*
 Copyright 2009 David Hall, Daniel Ramage

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/
package scalanlp.stage.source;

import java.io.File;

import scalanlp.pipes.Pipes;
import scalanlp.collection.LazyIterable;
import scalanlp.stage.{Parcel,Item,History};
import scalanlp.serialization._;
import scalanlp.ra.Cell;

/**
 * A TSV file acts as a source of Array[String] is a simple tab-delimited
 * file with no escaping.  The incoming strings are simply split on tab and
 * outgoing strings have their tabs replaced by a single space character.
 * 
 * @author dramage
 */
class TSVFile(path : String) extends File(path)
with LazyIterable[Seq[String]] with FileStreams {

  def write[V:TableWritable](table : V) = {
    val ps = printer();
    try {
      TSVTableSerialization.write(ps, table);
    } finally {
      ps.close();
    }
  }

  def read[V:TableReadable] : V = {
    val r = reader();
    try {
      TSVTableSerialization.read[V](r);
    } finally {
      reader.close;
    }
  }

  override def iterator =
    TSVTableSerialization.read[Iterator[List[String]]](this);

  override def toString =
    "TSVFile(\""+path+"\")";
}

/**
 * Companion object for TSVFiles that includes an implicit conversion
 * to a parcel and a format method that correctly escapes and joins
 * a sequence of strings.
 *
 * @author dramage
 */
object TSVFile {
  /** Named file in the current folder. */
  def apply(name : String)(implicit pipes : Pipes = Pipes.global) =
    new TSVFile(pipes.file(name).getPath);

  /** From file. */
  def apply(file : File) =
    new TSVFile(file.getPath);

  /** From file in directory. */
  def apply(file : File, name : String) =
    new TSVFile(new File(file, name).getPath);

  /** Calls file.asParcel. */
  implicit def TSVFileAsParcel(file : TSVFile) =
    Parcel(history = History.Origin(file.toString),
           meta = scalanlp.collection.immutable.DHMap() + (file : File),
           data = file.zipWithIndex.map(tup => Item(tup._2, tup._1)));

  /** Converts each whitespace character to a space. */
  def format(string : String) : String =
    string.replaceAll("\\s", " ");

  /** Formats the given sequence of strings as well-formed line of TSV. */
  def format(seq : Iterator[String]) : String =
    seq.map(format).mkString("\t");
}
