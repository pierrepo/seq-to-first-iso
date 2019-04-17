
"""Document containing the gestion of routes for the Flask application.

Uses flask_wtf for forms and seq-to-first-iso for computations.
"""

from pathlib import Path
import random
from threading import Thread

from flask import (flash,
                   redirect,
                   render_template,
                   request,
                   send_file,
                   session,
                   url_for)
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired
from wtforms import SelectMultipleField, SubmitField
from wtforms.widgets import CheckboxInput, ListWidget
from werkzeug.utils import secure_filename

from flask_app import app
from seq_to_first_iso import sequence_parser, seq_to_tsv


DICT_AMINO_ACIDS = {"A": "Alanine",
                    "C": "Cysteine",
                    "D": "Aspartic acid",
                    "E": "Glutamic acid",
                    "F": "Phenylalanine",
                    "G": "Glycine",
                    "H": "Histidine",
                    "I": "Isoleucine",
                    "K": "Lysine",
                    "L": "Leucine",
                    "M": "Methionine",
                    "N": "Asparagine",
                    "P": "Proline",
                    "Q": "Glutamine",
                    "R": "Arginine",
                    "S": "Serine",
                    "T": "Threonine",
                    "V": "Valine",
                    "W": "Tryptophan",
                    "Y": "Tyrosine",
                    }

# Use/create a folder named uploads in the current working directory.
UPLOAD_FOLDER = Path().cwd().joinpath("uploads")

if not UPLOAD_FOLDER.exists():
    UPLOAD_FOLDER.mkdir()

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
# Make key more secure ! + To remove if distributing the code ?
app.secret_key = "secret"
# Session to keep some informations.
app.config['SESSION_TYPE'] = 'filesystem'

threads = {}


class ThreadWithReturnValue(Thread):
    """Thread that returns a value.

    From : 'https://stackoverflow.com/questions/6893968/'.
    To get the value, once the thread is finished, call join().
    """

    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}, Verbose=None):
        """Initialize the thread."""
        Thread.__init__(self, group, target, name, args, kwargs)
        self._return = None

    def run(self):
        """Run the thread."""
        if self._target is not None:
            self._return = self._target(*self._args,
                                        **self._kwargs)

    def join(self, *args):
        """Get the value from the thread once finished."""
        Thread.join(self, *args)
        return self._return


class MultiCheckboxField(SelectMultipleField):
    """HTML field from flask_wtf with multiple checkboxes."""

    widget = ListWidget(prefix_label=False)
    option_widget = CheckboxInput()


class SequenceForm(FlaskForm):
    """Custom form for seq-to-first-iso."""

    # Add allowed file extensions ?
    upload = FileField("Sequences", validators=[FileRequired()])

    amino_acids = MultiCheckboxField("Amino acids",
                                     choices=[aa for
                                              aa in DICT_AMINO_ACIDS.items()]
                                     )

    submit = SubmitField('Submit')


@app.errorhandler(404)
def page_not_found(error):
    """Handle 404 errors to provide redirection."""
    return (render_template("error.html",
                            error=error,
                            info=f"URL: {request.path} "
                                 + "does not exist on the server"),
            404)


@app.route('/upload', methods=["GET", "POST"])
@app.route('/home', methods=["GET", "POST"])
@app.route('/index', methods=["GET", "POST"])
@app.route('/', methods=["GET", "POST"])
def wtfupload():
    """Route with WTForm to upload a file.

    Parse and verify the file format
    Creates a session with unlabelled amino acids
    Runs a thread from seq_to_tsv with random id
    """
    global threads
    form = SequenceForm()

    if form.validate_on_submit():
        # Make more verification on files ?

        # Save the file in UPLOAD_FOLDER.
        f = form.upload.data
        filename = secure_filename(f.filename)
        filepath = Path(app.config["UPLOAD_FOLDER"]).joinpath(filename)
        f.save(str(filepath))

        # Checked amino acids in the form.
        checked = [aa.data for aa in form.amino_acids if aa.checked]

        # Verification of file format.
        try:
            sequences, ignored_lines = sequence_parser(filepath)
        except UnicodeDecodeError as ude:
            return render_template("error.html",
                                   error=f"{ude}<br>"
                                         + "The file provided cannot be read",
                                   info="Make sure your file "
                                        + "is not a binary file")
        except Exception:
            return render_template("error.html",
                                   error="The file cannot be parsed")

        # The file can be read but no lines can be interpreted as sequences.
        if not sequences:
            return render_template("error.html",
                                   error="The file format is incorrect",
                                   info="No lines are only sequences")

        # Make it a var in session ?
        if ignored_lines:
            flash(f"{ignored_lines} lines ignored out of {len(sequences)}")

        session["amino_acids"] = ", ".join(checked)
        session["status"] = "Processing sequences"

        # Make sure that there can't be clashes after ?
        thread_id = random.randint(0, 10000)

        # Absolute path to the file (minus the .tsv).
        output_file = filepath.parent.joinpath(f"{filepath.stem}"
                                               + f"_{thread_id}")

        # Create a thread for seq_to_tsv and start it.
        threads[thread_id] = ThreadWithReturnValue(target=seq_to_tsv,
                                                   kwargs={"sequences":
                                                           sequences,
                                                           "output_file":
                                                           output_file,
                                                           "unlabelled_aa":
                                                           checked
                                                           })
        threads[thread_id].start()

        return redirect(url_for("progress", thread_id=thread_id))

    return render_template("upload.html", form=form)


@app.route("/progress/<int:thread_id>")
def progress(thread_id):
    """Progress page for a threaded task with the corresponding id."""
    global threads

    unlabelled_aa = session.get("amino_acids", [])
    status = session.get("status", "status unavailable")

    try:
        thread_running = threads[thread_id].isAlive()
    except KeyError as ke:
        return render_template("error.html",
                               error=f"Key error {ke}<br>"
                                     + "The process can't be accessed")

    if thread_running:
        return render_template("progress.html",
                               title="Loading",
                               unlabelled_aa=unlabelled_aa,
                               status=status)
    else:
        # Page with the link to download the file.
        return render_template("progress.html",
                               title="Done",
                               download_ready=True,
                               thread_id=thread_id)


@app.route("/results/<int:thread_id>")
def results(thread_id):
    """Page to return a file from a thread once computations are complete."""
    global threads

    try:
        output = threads[thread_id].join()
    except KeyError as ke:
        err_msg = f"Key error {ke}<br>The file can't be accessed<br>"
        info_msg = ("This might be due to the file being deleted "
                    + "or the server restarting<br>"
                    + "You need to launch another process again")
        return render_template("error.html",
                               error=err_msg,
                               info=info_msg)

    #output_name = output.stem + ".tsv"
    output_name = f"output_{thread_id}.tsv"
    output.to_csv(output_name, sep="\t", index=False)
    return send_file(filename_or_fp=output_name,
                     mimetype="text/tab-separated-values",
                     attachment_filename=output_name,
                     as_attachment=True)


if __name__ == '__main__':
    print("routes opened")
