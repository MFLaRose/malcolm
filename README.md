# malcolm
My pulsar data processing code for anyone to use and improve upon!

This code exists in two parts. The first is the actual code that I used while doing research, which was developed over many months of trial and error. The second is a set of jupyter notebooks and example files to guide a new user through the functionality of my code. 

NOTEBOOKS:
These are aimed at an auidience that is new to pulsar research and python. The difficulties of learning a programming language for the first time should not get in the way of good science, and it is to this end that I focussed myself on. It may be worth reading and executing a notebook multiple times to fully understand what is happening before diving into the scripts, especially if this is new for you. I also tried to explicity state when I was doing something because I didn't know a better method. 

There are also a few places where I have encouraged you to try editing the notebook yourself. This isn't entirely necessary, but will definitely aid in your understanding! 


NOTES ON RUNNING SCRIPTS AND NOTEBOOKS:
I have not included functionality for running either the scripts or notebooks multiple times. This is because we may want to make changes in the initial data processing steps, but may also want to be able to compare results later. The options were to either delete the old files, or check if they exist and save the new ones with a new naming scheme. I chose to do neither, as the first prevents us from making valuable comparisons, and the second will bury you in files. If you want to run these twice for one observation, I recommend moving the initial results to a new directory located above the code!

In order to start a notebook over from scratch, delete any files it has created and simply start from the beginning! :)

CRASH COURSE ON SCINTILLATION:

Welcome to the wonderful world of Interstellar Scintillation! This is by no means a simple subject, and it is my goal to prepare you for studying scintillation! After this section you will find a list of resources related to pulsars and scintillation. It is by no means exhasutive, but it should be a useful reference. I will try to describe scintillation in as many ways as I understand it, starting from simple metaphors, and then talking about a little physics. 

If you've ever looked up at the stars at night you will have certainly noticed that stars twinkle, but planets are a consistent source of light. When you stare at a twinkling star for a while, you'll notice it get sporadically brighter and dimmer, and occasionally even shift a little in color. What if I told you that you could measure properties about the atmosphere with this twinkling? This is essentially what we are doing when we study scintillation! As radio light from the pulsar travels through the galaxy, it passes through the interstellar medium, which is largely composed of ions. This diffuse gas does to a pulsar's signal what the atmosphere does to a star's light! So when we are observing a pulsar, scintillation shows up as a characteristic fluctuation in intensity over frequency and time. This is why this is a useful analogy! As we 'stare' at a pulsar, we are essentially watching it get sporadically brighter and dimmer, and occasionally even shift a little in frequency (color!). 

Now it's time for some physics! A simple model of the pulsar-ISM-earth system treats the pulsar as a point source whose signal passes through a 'screen' of gas. As we know from elementary physics, light passes through different media with different speeds. If the interstellar gas were uniform (or if you want to sound fancy, homogenous), and that is to say that the index of refraction is constant throughout the medium, we would see the pulsar's signal systematically shifted in space, much like the classic pencil in a glass of water illusion. However, we live in a dynamic galaxy filled with stars and rocks and gas flying in all sorts of directions. Essentially, the gas in space is turbulent, which means it is not a uniform medium! So instead of the signal being shifted, we see a pattern of light and dark spots where the light waves constructively and destructively interfered, respectively. This pattern is what we are looking at when we create a dynamic spectrum! To measure properties about the scintillation itself, we measure the size of these blobs (also known as scintles) in frequency and time, which we call the scintillation bandwidth and scintillation timescale. 

And now to make things a little more complicated. As I said before, we live in a dynamic galaxy. Not only is the gas in the ISM turbulent, it's also generally moving in the same direction. And because that wasn't enough, pulsars are also moving through space, and so is the earth! And so the pattern of light and dark spots also moves across the observer (us!) systematically, at a rate determined by the relative motions of the three parts of the system. I like to think of this problem as a problem with three variables. For example, if we know something about the relative velocities of the pulsar and earth, we can use observations of scintillation to say something about the motion and structure of the ISM. 

So there are really two reasons to study scintillation. The first is to learn more about the structure of our galaxy. We can do this by establishing a model for the scattering medium in the galaxy (in this case, we use the NE2001 model of galactic electron density) and making observations of scintillation. Differences in the model's predictions for scintillation paramters and the observed paramters indicate places where the model can (and should) be improved upon. There are even some papers that discuss using scintillometry (the measure of scintillation) to probe the ISM on *solar system sized scales*! So making really precise measurements of scintillation helps us make really precise models of the galaxy!
The second reason is to learn more about pulsars themselves. As I said before, I like to think about this like a problem with three variables. If we have knowledge about two of them, we can solve for the third. There are many areas of pulsar science that would benefit from more accurate measurements of their motion through space. An active area of study that comes to mind is pulsar formation. Pulsars form at the end of a massive star's life during a supernova, where the failure of nuclear fusion causes massive compression of a star's material. This supernova explosion can impart a 'kick' velocity on the pulsar, much like debris flying away from an explosion. So if we can constrain the kick velocity of a pulsar, we can say something about the energetic event that formed it!

Thank you for coming to my TED talk! If you want to read more about pulsars and scintillation (and you really should, this is cool stuff!), please look through the links below at your own liesure! 

USEFUL LINKS:

Primers:
* https://en.wikipedia.org/wiki/Pulsar
* https://en.wikipedia.org/wiki/Neutron_star
* https://en.wikipedia.org/wiki/Interstellar_medium
* https://en.wikipedia.org/wiki/Diffraction
* http://pulsarsearchcollaboratory.com/online-workshop/

Important Info:
* http://www.apc.univ-paris7.fr/Downloads/astrog/djannati/Flower/HandbookOfPulsarAstronomy_Lorimer-Kramer.pdf
* https://iopscience.iop.org/article/10.1086/306358/fulltext/38012.text.html
* http://articles.adsabs.harvard.edu/pdf/1986ApJ...311..183C
* https://arxiv.org/pdf/astro-ph/0207156.pdf

Deeper Dive:
* https://iopscience.iop.org/article/10.1086/319133/fulltext/005821.text.html
* https://arxiv.org/pdf/1601.04490.pdf
* https://arxiv.org/pdf/1302.1897.pdf
* https://www.cita.utoronto.ca/~fkirsten/lib/exe/fetch.php?media=start:workshop:j-p-macquart-bonnscintorthodoxytalk.pdf
* https://academic.oup.com/mnras/article/354/1/43/963584
* https://iopscience.iop.org/article/10.1086/307181/pdf





